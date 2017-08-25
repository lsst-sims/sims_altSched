from __future__ import division

from graphics.GraphicalMonitor import GraphicalMonitor
#from graphics.GraphicalMonitor3D import GraphicalMonitor3D

import numpy as np
from minis.Scheduler import Scheduler
from Visit import PROP_WFD, PROP_DD
import Config
from lsst.sims.speedObservatory import sky
from lsst.sims.speedObservatory import Telescope
from SkyMap import SkyMap
import time
import sys
import os
from multiprocessing import Pool
from datetime import datetime

from matplotlib import pyplot as plt
from SummaryPlots import SummaryPlots

from lsst.sims.ocs.kernel.time_handler import TimeHandler
from lsst.sims.ocs.environment.cloud_model import CloudModel

from lsst.sims.speedObservatory.utils import unix2mjd

trackMap = True
showDisp = True and trackMap
graphicsMode = "filters"
saveMovie = False
showSummaryPlots = True and trackMap
clearDisplayNightly = True
writeCsv = False

resScale = 3
runName = "9011_4vo_maxf_lax_3set_97"

assert(not (showDisp and not trackMap))

class Simulator:
    def __init__(self):
        pass

    def time(self):
        return self.curTime

    def run(self, tel):
        self.tel = tel

        # read in downtime nights TODO use SOCS/speedObservatory for this
        downtimeNights = self.parseDowntime()

        # initialize the scheduler
        self.sched = Scheduler(telescope=self.tel, context=self)

        # make a TimeHandler to give to the CloudModel
        # TODO use sims_speedObservatory instead
        dateFormat = "%Y-%m-%d"
        startDatetime = datetime.utcfromtimestamp(Config.surveyStartTime)
        timeHandler = TimeHandler(datetime.strftime(startDatetime, dateFormat))
        self.cloudModel = CloudModel(timeHandler)
        # load the cloud database
        self.cloudModel.initialize()

        # resScale is proportional to the resolution (and therefore the speed)
        if trackMap:
            self.skyMap = SkyMap(telescope=self.tel, resScale=resScale)
        if showDisp:
            self.display = GraphicalMonitor(skyMap=self.skyMap, mode=graphicsMode)

        # track slew times for simplest summary statistics
        self.slewTimes = []

        # write the header to the output file if writeCsv flag is set
        if writeCsv:
            self.outFile = open("results/" + runName + "/" + runName + ".csv", "w")
            self.outFile.write("time,prop,ra,dec,filter\n")

        # run the survey!
        # these variables keep track of wasted time
        self.fieldsRisingWasteTime = 0
        self.earlyNightEndWasteTime = 0
        for nightNum in range(Config.surveyNumNights):
            print "Night:", nightNum, "\r",
            sys.stdout.flush()
            if nightNum in downtimeNights:
                continue

            self.simulateNight(nightNum)
            self.sched.notifyNightEnd()

        print "time wasted waiting for fields to rise:", self.fieldsRisingWasteTime
        print "time wasted when sched ran out of visits:", self.earlyNightEndWasteTime
        if writeCsv:
            self.outFile.close()

        # we're done with the simulation now
        # show summary plots if requested
        if showSummaryPlots:
            self._outputSummaryStats()

    def simulateNight(self, nightNum):
        twilStart = sky.twilStart(Config.surveyStartTime, nightNum)
        twilEnd   = sky.twilEnd(Config.surveyStartTime, nightNum)

        # start out the simulation at the beginning of the night
        self.curTime = twilStart

        # prevI is the last value of i when we updated the display
        prevI = 0

        # prevFilter is the filter of the last visit
        (prevFilter, prevAlt, prevAz) = ('', -1, -1)
        wasDomeJustClosed = False
        for i, visit in enumerate(self.sched.scheduleNight(nightNum)):
            perNight, deltaI = self.getUpdateRate(nightNum, i)

            # stop if the night is over
            if self.curTime >= twilEnd:
                break

            # skip forward in time if there are clouds
            deltaT = self.curTime - Config.surveyStartTime
            cloudCover = self.cloudModel.get_cloud(deltaT)
            timeBeforeDomeClose = self.curTime
            while cloudCover > Config.maxCloudCover and self.curTime <= twilEnd:
                self.curTime += 600
                deltaT = self.curTime - Config.surveyStartTime
                cloudCover = self.cloudModel.get_cloud(deltaT)

            # if the dome closed due to clouds, get a fresh visit from the
            # scheduler by continuing
            if self.curTime > timeBeforeDomeClose:
                # let sched know that the dome closed for a while
                self.sched.notifyDomeClosed(self.curTime - timeBeforeDomeClose)
                wasDomeJustClosed = True
                # continue to get a new visit that isn't stale
                continue

            # if visit is None, that means the scheduler has nowhere to point
            # at the moment
            if visit is None:
                self.curTime += 30
                self.fieldsRisingWasteTime += 30
            else:
                alt, az = sky.radec2altaz(visit.ra, visit.dec, self.curTime)
                # make sure this az is a valid place to look
                if alt < self.tel.minAlt or alt > self.tel.maxAlt:
                    print "invalid alt (", np.degrees(alt), "deg) night", nightNum
                    continue

                # figure out how far we have to slew
                if i > 0 and not wasDomeJustClosed:
                    slewTime = self.tel.calcSlewTime(prevAlt, prevAz, prevFilter,
                                                     alt, az, visit.filter,
                                                     laxDome = Config.laxDome)
                    self.curTime += slewTime
                    self.slewTimes.append(slewTime)

                # notify the skyMap of the visit
                # (the time of the visit is the time after the slew is over)
                if trackMap:
                    self.skyMap.addVisit(visit, self.curTime)

                if writeCsv:
                    assert(visit.prop == PROP_WFD or visit.prop == PROP_DD)
                    prop = "wfd" if visit.prop == PROP_WFD else "dd"
                    self.outFile.write(str(unix2mjd(self.curTime)) + "," +
                                       prop + "," +
                                       str(visit.ra) + "," +
                                       str(visit.dec) + "," +
                                       visit.filter + "\n")

                # add the exposure time of this visit to the current time
                self.curTime += visit.expTime
                self.curTime += Config.visitOverheadTime

                # let the scheduler know we "carried out" this visit
                self.sched.notifyVisitComplete(visit, self.curTime)

                (prevAlt, prevAz, prevFilter) = (alt, az, visit.filter)

            wasDomeJustClosed = False

            # now that we've added the visit (if there was one),
            # update the display
            if showDisp and not perNight and i - prevI >= deltaI:
                self.display.updateDisplay(self.skyMap, self.curTime)
                prevI = i
                # save each frame if the saveMovie flag is set
                if saveMovie:
                    self.display.saveFrame("images/pygame/%07d.png" % i)

        # the night is over

        # calculate how much time we waste at the end of the night
        self.earlyNightEndWasteTime += twilEnd - self.curTime
        # show the sky scroll by during the wasted time
        if showDisp and not perNight:
            while self.curTime < twilEnd:
                self.display.updateDisplay(self.skyMap, self.curTime)
                self.curTime += 30

        # update the display with the completed night if perNight is true
        if showDisp and perNight:
            self.display.updateDisplay(self.skyMap, self.curTime)
            if saveMovie:
                self.display.saveFrame("images/pygame/%07d.png" % i)

        # clear the skyMap so the next updateDisplay won't have tonight's visits
        if showDisp and clearDisplayNightly:
            self.skyMap.clear()

    def _outputSummaryStats(self):
        """ Displays some example summary statistics

        These have mostly been superceded by MAF, but could potentially
        be useful still as a quick assessment of a run.
        """

        avgRevisitTimes = self.skyMap.getAvgRevisitMap()
        plt.figure("revisit times")
        plt.title("Average Revisit Times (in minutes)")
        plt.imshow(avgRevisitTimes/60)
        plt.colorbar()

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

        plt.figure("% Visits not accompanied by a revisit within 45 minutes (real hrs)")
        percentLonelyMap = self.skyMap.getLonelinessMap(cutoffMins=45)
        plt.title("Fraction of visits with no revisit within 45 minutes")
        plt.imshow(percentLonelyMap)
        plt.clim(0,0.4)
        plt.colorbar()

        for percentile in [10, 50, 75, 90, 95, 99]:
            plt.figure(str(percentile) + "th percentile revisit time")
            percentileMap = self.skyMap.getPercentileMap(percentile)
            plt.title(str(percentile) + "th percentile revisit time (in days)")
            plt.imshow(percentileMap / 3600 / 24)
            plt.clim(0,7)
            plt.colorbar()

        plotter = SummaryPlots(self.skyMap, slewTimes = self.slewTimes)

        print "avg slew time", np.mean(self.slewTimes), "seconds"
        print "median slew time", np.median(self.slewTimes), "seconds"
        plotter.slewHist()
        
        sortedTimes = np.sort(self.slewTimes)
        cum = np.cumsum(sortedTimes)
        print "total cumulative slew time: ", cum[-1]
        print "rank @ half total cum / # slews", np.searchsorted(cum, cum[-1]/2) / len(cum)
        plotter.slewRankCum()
        
        plotter.revisitHist()
        
        plotter.dAirmassCum()
        plotter.dAirmassContour()
        plotter.zenithAngleContour()
        plotter.airmassHist()

        plotter.show()

    @staticmethod
    def getUpdateRate(nightNum, i):
        # returning (False, 1) makes the display just show every visit
        return (False, 1)

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

        return (perNight, deltaI)

    # TODO use sims_speedObservatory instead
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



def runDefaultSim():
    sim = Simulator()
    # use a telescope with default parameters
    tel = Telescope()
    if not os.path.isdir("results/" + runName):
        os.mkdir("results/" + runName)
    return sim.run(tel)

if __name__ == "__main__":
    runDefaultSim()
