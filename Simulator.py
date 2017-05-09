from __future__ import division

from graphics.GraphicalMonitor import GraphicalMonitor
#from graphics.GraphicalMonitor3D import GraphicalMonitor3D

import numpy as np
from minis.Scheduler import Scheduler
import Config
import AstronomicalSky
from Telescope import Telescope
from minis.SkyMap import SkyMap
import time
import sys
from multiprocessing import Pool

from matplotlib import pyplot as plt
from SummaryPlots import SummaryPlots

trackMap = True
showDisp = True
saveMovie = False
plotAzes = False
showSummaryPlots = True
clearDisplayNightly = False

# TODO plot azimuth vs time
# change filters every revisit set (revisits 1 hour, change every 2 hours)
# add in settle time

class Simulator:
    def __init__(self):
        self.startTime = Config.surveyStartTime
        self.curTime = self.startTime

    @staticmethod
    def getUpdateRate(i):
        # decide how often to update the display so we can get
        # a speeding up effect
        perNight = False
        deltaI = 0
        if 0 <= i < 1000:
            deltaI = 1
        if 1000 <= i < 1500:
            deltaI = int((i - 1000)/500*10) + 1
        if 1500 <= i < 10000:
            deltaI = 10
        if 10000 <= i < 50000:
            deltaI = int((i - 10000)/50000*900) + 10
        if 50000 <= i:
            perNight = True

        return (True, 1)
        return (perNight, deltaI)

    def run(self, tel):
        totalNVisits = 0
        sched = Scheduler(telescope=tel, context=self)

        # resScale is proportional to the resolution (and therefore the speed)
        if trackMap:
            skyMap = SkyMap(telescope=tel, resScale=2)
        if showDisp:
            display = GraphicalMonitor(skyMap=skyMap)
       
        nightNum = 0
        isNightYoung = True

        slewTimes = []
        revisitTimes = []

        # keep track of the alt and az of each visit
        alts = []
        azes = []
        obsDecs = []
        if plotAzes:
            plt.clf()
            plt.xlabel("visit number")
            plt.ylabel("azimuth (degrees)")
            plt.ylim(0, 360)
            plt.ion()
            plt.show()

        i = 0
        prevI = i
        for visit in sched.schedule():

            perNight, deltaI = self.getUpdateRate(i)
            if showDisp and ((    perNight and isNightYoung) or
                             (not perNight and i - prevI >= deltaI)):
                display.updateDisplay(skyMap, self.curTime)
                prevI = i
                if saveMovie:
                    display.saveFrame("images/pygame/%07d.png" % i)

            if isNightYoung:
                print "Night:", nightNum, "\r",
                sys.stdout.flush()
            # if visit is None, that means the scheduler ran out of places
            # to point for the night
            if visit is not None:
                totalNVisits += 1
                # keep track of the alt and az of each visit
                radec = np.array([[visit.ra, visit.dec]])
                obsDecs.append(visit.dec)
                altaz = AstronomicalSky.radec2altaz(radec, self.time())[0]
                alts.append(altaz[0])
                azes.append(altaz[1])
                if not isNightYoung:
                    slewTime = tel.calcSlewTime(prevAltaz, altaz)

                # notify the display of the visit
                if trackMap:
                    skyMap.addVisit(visit, self.curTime)

            # plot the azimuth at each time step
            if plotAzes:
                plt.scatter(range(len(azes)), np.degrees(azes))
                plt.draw()
                plt.pause(0.01)

            if visit is None:
                self.curTime += 30
            else:
                # add the exposure time of this visit to the current time
                expTime = AstronomicalSky.getExpTime(visit.ra, visit.dec)
                self.curTime += expTime
                if not isNightYoung:
                    self.curTime += slewTime
                    slewTimes.append(slewTime)
                sched.notifyVisitComplete(visit, self.time())

                # keep track of revisit times
                visitPair = visit.visitPair
                if visitPair.visit1.isComplete and visitPair.visit2.isComplete:
                    revisitTimes.append(np.abs(visitPair.visit1.timeOfCompletion - \
                                               visitPair.visit2.timeOfCompletion))
                prevAltaz = altaz

            if isNightYoung and showDisp and clearDisplayNightly:
                skyMap.clear()
            # process the end of the night if necessary
            if self.curTime < AstronomicalSky.nightEnd(nightNum):
                isNightYoung = False
            else:
                self.curTime = AstronomicalSky.nightStart(nightNum + 1)
                nightNum += 1
                isNightYoung = True

            i += 1
            if nightNum > (365/12):
                break

        if not showSummaryPlots:
            return

        avgRevisitTimes = skyMap.getAvgRevisitMap()
        plt.figure("revisit times")
        plt.title("Average Revisit Times (in minutes)")
        plt.imshow(avgRevisitTimes/60)
        plt.colorbar()

        revisitTimesMap = skyMap.getRevisitMap()
        allRevisitTimes = []
        for pix in revisitTimesMap.flatten():
            # some entries might be zeros
            if isinstance(pix, list):
                allRevisitTimes += pix
        plt.figure("Revisit Time Histogram")
        plt.hist(np.array(allRevisitTimes)/3600, 300, cumulative=True, normed=True, range=[0,2])
        plt.title("Per-Pixel Revisit Times")
        plt.xlabel("Time (hours)")
        plt.ylabel("Cumulative number of visits (normalized)")

        plt.figure("% Visits not accompanied by a revisit within 45 minutes")
        percentLonelyMap = skyMap.getLonelinessMap(cutoffMins=45)
        plt.title("Percent of visits with no revisit within 45 minutes")
        plt.imshow(percentLonelyMap)
        plt.clim(0,0.4)
        plt.colorbar()

        for percentile in [10, 50, 75, 90, 95, 99]:
            plt.figure(str(percentile) + "th percentile revisit time")
            percentileMap = skyMap.getPercentileMap(percentile)
            plt.title(str(percentile) + "th percentile revisit time (in days)")
            plt.imshow(percentileMap / 3600 / 24)
            plt.clim(0,7)
            plt.colorbar()

        azes = np.array(azes) % (2*np.pi)

        plotter = SummaryPlots(alts = alts, azes = azes, decs = obsDecs,
                               slewTimes = slewTimes, revisitTimes = revisitTimes)

        print "avg slew time", np.mean(slewTimes), "seconds"
        print "median slew time", np.median(slewTimes), "seconds"
        plotter.slewHist()
        
        sortedTimes = np.sort(slewTimes)
        cum = np.cumsum(sortedTimes)
        print "total cumulative slew time: ", cum[-1]
        print "rank @ half total cum / # slews", np.searchsorted(cum, cum[-1]/2) / len(cum)
        plotter.slewRankCum()
        
        print "avg revisit time", np.mean(np.array(revisitTimes)/60), "minutes"
        plotter.revisitHist()
        
        plotter.dAirmassCum()
        plotter.dAirmassContour()
        plotter.airmassHist()

        plotter.show()

    def time(self):
        return self.curTime

def runSimWorker(tel):
    sim = Simulator()
    return sim.run(tel)

def runDefaultSim():
    sim = Simulator()
    # use a telescope with default parameters
    tel = Telescope()
    return sim.run(tel)

if __name__ == "__main__":
    runDefaultSim()
    exit()
    domAzMaxSpeed = 1.5
    domAzMaxSpeeds = np.arange(1.5, 2.1, 0.1)
    domAltMaxSpeeds = np.arange(1.75, 4, 0.25)
    domAltMaxSpeed = 1.75
    #domAltMaxSpeeds = np.arange(1.75, 2.25, 0.25)
    settleTimes = np.arange(2.8,3.1,0.1)
    #settleTimes = np.arange(2.5,3,0.25)
    results = np.zeros((len(domAzMaxSpeeds), len(settleTimes)))
    print
    print "tot:", len(domAzMaxSpeeds)
    print
    for i, domAzMaxSpeed in enumerate(domAzMaxSpeeds):
        print
        print "i:", i
        print
        p = Pool(8)
        tels = [Telescope() for settleTime in settleTimes]
        for tel, settleTime in zip(tels, settleTimes):
            tel.domAltMaxSpeed = domAltMaxSpeed
            tel.domAzMaxSpeed = domAzMaxSpeed
            tel.settleTime = settleTime
        row = p.map(runSimWorker, tels)
        results[-i-1,:] = row
        #for j, settleTime in enumerate(settleTimes):
        #    sim = Simulator()
        #    results[-i-1, j] = sim.run(domAltMaxSpeed, domAzMaxSpeed, settleTime)
    print results.tolist()
    plt.matshow(results, extent=[settleTimes.min(), settleTimes.max(), domAzMaxSpeeds.min(), domAzMaxSpeeds.max()])
    ax = plt.gca()
    ax.set_xticks(settleTimes[::5])
    ax.set_yticks(domAzMaxSpeeds[::5])
    plt.colorbar()
    plt.title("Total Number of Visits (six-month sims)")
    plt.xlabel("Settle Time (secs)")
    plt.ylabel("Dome Azimuth Max Speed (deg/sec)")
    plt.show()
