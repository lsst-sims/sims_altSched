from __future__ import division

from graphics.GraphicalMonitor import GraphicalMonitor
#from graphics.GraphicalMonitor3D import GraphicalMonitor3D

import numpy as np
from minis.Scheduler import Scheduler
import Config
import AstronomicalSky
import Telescope
from minis.SkyMap import SkyMap
import time
import sys

from matplotlib import pyplot as plt
from SummaryPlots import SummaryPlots

showDisp = False
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

        return (False, 1)
        return (perNight, deltaI)

    def run(self):
        sched = Scheduler(context=self)

        # resScale is proportional to the resolution (and therefore the speed)
        if showDisp:
            skyMap = SkyMap(resScale=6)
            display = GraphicalMonitor(context=self, skyMap=skyMap)
       
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
            if isNightYoung:
                print "Night:", nightNum, "\r",
                sys.stdout.flush()
            # if visit is None, that means the scheduler ran out of places
            # to point for the night
            if visit is not None:
                # keep track of the alt and az of each visit
                radec = np.array([[visit.ra, visit.dec]])
                obsDecs.append(visit.dec)
                altaz = AstronomicalSky.radec2altaz(radec, self.time())[0]
                alts.append(altaz[0])
                azes.append(altaz[1])
                if not isNightYoung:
                    slewTime = Telescope.calcSlewTime(prevAltaz, altaz)

                # notify the display of the visit
                if showDisp:
                    skyMap.addVisit(visit, self.curTime)

            # plot the azimuth at each time step
            if plotAzes:
                plt.scatter(range(len(azes)), np.degrees(azes))
                plt.draw()
                plt.pause(0.01)

            perNight, deltaI = self.getUpdateRate(i)

            if showDisp and ((    perNight and isNightYoung) or
                             (not perNight and i - prevI >= deltaI)):
                display.updateDisplay(skyMap)
                prevI = i
                if saveMovie:
                    display.saveFrame("images/pygame/%07d.png" % i)

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

            # process the end of the night if necessary
            if self.curTime < AstronomicalSky.nightEnd(nightNum):
                isNightYoung = False
            else:
                self.curTime = AstronomicalSky.nightStart(nightNum + 1)
                nightNum += 1
                isNightYoung = True
                if showDisp and clearDisplayNightly:
                    skyMap.clear()

            i += 1
            if i > 20000:
                print
                break


        if not showSummaryPlots:
            return

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

if __name__ == "__main__":
    sim = Simulator()
    sim.run()
