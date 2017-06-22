from __future__ import division

from graphics.GraphicalMonitor import GraphicalMonitor
#from graphics.GraphicalMonitor3D import GraphicalMonitor3D

import numpy as np
from minis.Scheduler import Scheduler
import Config
from lsst.sims.speedObservatory import sky
from lsst.sims.speedObservatory import Telescope
from SkyMap import SkyMap
import time
import sys
from multiprocessing import Pool
from datetime import datetime

from matplotlib import pyplot as plt
from SummaryPlots import SummaryPlots

from lsst.sims.ocs.kernel.time_handler import TimeHandler
from lsst.sims.ocs.environment.cloud_model import CloudModel

trackMap = True
showDisp = True and trackMap
saveMovie = False
showSummaryPlots = True and trackMap
clearDisplayNightly = True
writeCsv = False

assert(not (showDisp and not trackMap))

surveyYears = 10

class Simulator:
    def __init__(self):
        self.curTime = sky.nightStart(Config.surveyStartTime, 1)

    @staticmethod
    def getUpdateRate(nightNum, i):
        # TODO hacky -- assume 1000 visits/night to create an "effective" i
        i = nightNum * 1000 + i
        # decide how often to update the display so we can get
        # a speeding up effect
        perNight = False
        deltaI = 0
        speed1 = 10
        speed2 = 900
        cut1 = 1000
        cut2 = 1500
        cut3 = 10000
        cut4 = 50000
        if 0 <= i < cut1:
            deltaI = 1
        if cut1 <= i < cut2:
            deltaI = int((i - cut1)/(cut2-cut1)*speed1) + 1
        if cut2 <= i < cut3:
            deltaI = speed1
        if cut3 <= i < cut4:
            deltaI = int((i - cut3)/(cut4-cut3)*speed2) + speed1
        if cut4 <= i:
            perNight = True

        return (False, 1)
        return (perNight, deltaI)

    def parseDowntime(self, schedFileName="schedDown.conf",
                            unschedFileName="unschedDown.conf"):
        def parseSchedFile(filename):
            downtimeNights = set([])
            with open(filename) as schedFile:
                for line in schedFile:
                    line = line.strip()
                    if len(line) == 0 or line[0] == "#":
                        continue
                    line = line.split("#")[0].strip()
                    param, value = map(str.strip, line.split("="))
                    if param == "startNight":
                        startNight = int(value)
                    if param == "duration":
                        # relying on startNight having already been set
                        duration = int(value)
                        downtimeNights.update(range(startNight, startNight + duration))
            return downtimeNights

        schedDownNights = parseSchedFile(schedFileName)
        unschedDownNights = parseSchedFile(unschedFileName)
        return schedDownNights | unschedDownNights

    def run(self, tel):
        self.tel = tel

        # read in downtime nights TODO use SOCS/speedObservatory for this
        downtimeNights = self.parseDowntime()

        # initialize the scheduler
        self.sched = Scheduler(telescope=self.tel, context=self)

        # make a TimeHandler to give to the CloudModel
        dateFormat = "%Y-%m-%d"
        startDatetime = datetime.utcfromtimestamp(Config.surveyStartTime)
        timeHandler = TimeHandler(datetime.strftime(startDatetime, dateFormat))
        self.cloudModel = CloudModel(timeHandler)
        # load the cloud database
        self.cloudModel.initialize()

        # resScale is proportional to the resolution (and therefore the speed)
        if trackMap:
            self.skyMap = SkyMap(telescope=self.tel, resScale=6)
        if showDisp:
            self.display = GraphicalMonitor(skyMap=self.skyMap, mode="nvisits")
       
        self.slewTimes = []
        self.revisitTimes = []

        # write the header to the output file if writeCsv flag is set
        if writeCsv:
            outFile = open("results/downtime_nightlen.csv", "w")
            outFile.write("time (unix), ra (rads), dec (rads), filter\n")

        # run the survey!
        numSimulatedNights = int(surveyYears * 365.25)
        for nightNum in range(numSimulatedNights):
            print "Night:", nightNum, "\r",
            sys.stdout.flush()
            if nightNum in downtimeNights:
                continue
            self.simulateNight(nightNum)

            # the clearDisplayNightly flag indicates we should clear
            # the skymap at the end of each night
            if showDisp and clearDisplayNightly:
                self.skyMap.clear()

        if writeCsv:
            outFile.close()

        # we're done with the simulation now unless we have to show
        # the summary plots
        if not showSummaryPlots:
            return

        """
        avgRevisitTimes = self.skyMap.getAvgRevisitMap()
        plt.figure("revisit times")
        plt.title("Average Revisit Times (in minutes)")
        plt.imshow(avgRevisitTimes/60)
        plt.colorbar()
        """

        """
        revisitTimesMap = self.skyMap.getRevisitMap()
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
        """

        plt.figure("% Visits not accompanied by a revisit within 45 minutes (real hrs)")
        percentLonelyMap = self.skyMap.getLonelinessMap(cutoffMins=45)
        plt.title("Fraction of visits with no revisit within 45 minutes")
        plt.imshow(percentLonelyMap)
        plt.clim(0,0.4)
        plt.colorbar()

        """
        for percentile in [10, 50, 75, 90, 95, 99]:
            plt.figure(str(percentile) + "th percentile revisit time")
            percentileMap = self.skyMap.getPercentileMap(percentile)
            plt.title(str(percentile) + "th percentile revisit time (in days)")
            plt.imshow(percentileMap / 3600 / 24)
            plt.clim(0,7)
            plt.colorbar()
        """


        plotter = SummaryPlots(self.skyMap, slewTimes = self.slewTimes)

        print "avg slew time", np.mean(self.slewTimes), "seconds"
        print "median slew time", np.median(self.slewTimes), "seconds"
        #plotter.slewHist()
        
        #sortedTimes = np.sort(self.slewTimes)
        #cum = np.cumsum(sortedTimes)
        #print "total cumulative slew time: ", cum[-1]
        #print "rank @ half total cum / # slews", np.searchsorted(cum, cum[-1]/2) / len(cum)
        #plotter.slewRankCum()
        
        print "avg revisit time", np.mean(np.array(self.revisitTimes)/60), "minutes"
        #plotter.revisitHist()
        
        plotter.dAirmassCum()
        #plotter.dAirmassContour()
        plotter.zenithAngleContour()
        #plotter.airmassHist()

        plotter.show()

    def simulateNight(self, nightNum):
        nightStart = sky.nightStart(Config.surveyStartTime, nightNum)
        nightEnd   = sky.nightEnd(Config.surveyStartTime, nightNum)

        # start out the simulation at the beginning of the night
        self.curTime = nightStart

        # prevI is the last value of i when we updated the display
        prevI = 0

        # prevFilter is the filter of the last visit
        prevFilter = ''
        for i, visit in enumerate(self.sched.scheduleNight(nightNum)):
            perNight, deltaI = self.getUpdateRate(nightNum, i)

            # if visit is None, that means the scheduler ran out of places
            # to point for the night
            if visit is None:
                self.curTime += 30
            else:
                alt, az = sky.radec2altaz(visit.ra, visit.dec, self.time())
                # make sure this az is a valid place to look
                if alt < 0:
                    raise RuntimeError("Can't look at the ground! " +
                                       "visit=" + visit)
                if alt > self.tel.maxAlt:
                    #print "Warning: tried to observe in zenith avoid zone"
                    continue

                # figure out how far we have to slew
                if i > 0:
                    slewTime = self.tel.calcSlewTime(prevAlt, prevAz, prevFilter,
                                                     alt, az, visit.filter)

                # notify the skyMap of the visit
                if trackMap:
                    self.skyMap.addVisit(visit, self.curTime)

                # add the exposure time of this visit to the current time
                expTime = sky.getExpTime(visit.ra, visit.dec)
                self.curTime += expTime
                if i > 0:
                    self.curTime += slewTime
                    self.slewTimes.append(slewTime)
                if writeCsv:
                    outFile.write(str(self.time()) + "," + \
                                  str(visit.ra) + "," + \
                                  str(visit.dec) + "," + \
                                  visit.filter + "\n")
                # let the scheduler know we "carried out" this visit
                self.sched.notifyVisitComplete(visit, self.time())

                # keep track of revisit times
                visit1 = visit.visitPair.visit1
                visit2 = visit.visitPair.visit2
                if visit1.isComplete and visit2.isComplete:
                    self.revisitTimes.append(np.abs(visit1.timeOfCompletion - \
                                                    visit2.timeOfCompletion))
                prevAlt = alt
                prevAz = az
                prevFilter = visit.filter


            # now that we've added the visit (if there was one),
            # update the display
            if showDisp and ((    perNight and i == 0) or
                             (not perNight and i - prevI >= deltaI)):
                self.display.updateDisplay(self.skyMap, self.curTime)
                prevI = i
                # save each frame if the saveMovie flag is set
                if saveMovie:
                    self.display.saveFrame("images/pygame/%07d.png" % i)

            # skip forward in time if there are clouds
            deltaT = self.curTime - Config.surveyStartTime
            cloudCover = self.cloudModel.get_cloud(deltaT)
            while cloudCover > Config.maxCloudCover and self.curTime <= nightEnd:
                self.curTime += 600
                self.display.updateDisplay(self.skyMap, self.curTime)
                deltaT = self.curTime - Config.surveyStartTime
                cloudCover = self.cloudModel.get_cloud(deltaT)
            # process the end of the night if necessary
            if self.curTime > nightEnd:
                break
        self.sched.notifyNightEnd()

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
