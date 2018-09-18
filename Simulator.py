from __future__ import division
from __future__ import print_function

from graphics.GraphicalMonitor import GraphicalMonitor

import numpy as np
from minis.Scheduler import Scheduler
from minis.XTwilScheduler import XTwilScheduler
from Visit import PROP_WFD, PROP_DD
import config
from lsst.sims.speedObservatory import sky
from lsst.sims.speedObservatory import Telescope
from SkyMap import SkyMap
import sys
import os
from datetime import datetime
from config import NORTH, SOUTH
import time

from matplotlib import pyplot as plt
from SummaryPlots import SummaryPlots

from lsst.sims.ocs.kernel.time_handler import TimeHandler
# TODO
from lsst.sims.ocs.environment.cloud_model import CloudModel

from lsst.sims.speedObservatory.utils import unix2mjd

# flag to control whether a SkyMap instance is used to track how many times
# each sky pixel is observed
trackMap = False

# flag to control whether a window pops up to visualize the progress of the
# scheduler
showDisp = True and trackMap

# flag to control whether we're simulating extended twilight time
doXTwil = True

# which mode to run the graphics in (it's passed as the argument to the
# GraphicalMonitor constructor). Options are "nvisits" -- visualize the total
# number of visits to each pixel -- and "filters" -- visualize the last
# filter used to observe each pixel
graphicsMode = "filters"

# flag to control whether to save each frame of the visualization to disk
saveMovie = False and showDisp

# location on disk to save each frame
movieDir = "images/pygame"

# flag to control whether summary plots are displayed after the simulation
# is complete (whether or not to call _outputSummaryStats())
showSummaryPlots = True and trackMap

# flag to control whether or not to clear the visualization after each night
clearDisplayNightly = True

# flag to control whether the observation sequence is written to disk
writeCsv = False

# the resolution of the visualization (using 2x will give 4 times as many pixels
# as using x)
resScale = 4

# the name of the run. The csv output is saved to
# ./results/`runName`/`runName`.csv
runName = "9013_1vo_maxf_nolax_3set_98"

class Simulator:
    def time(self):
        """ Get the current simulated time

        Returns
        -------
        The current simulated time as a unix timestamp
        """
        return self.curTime

    def run(self, tel):
        """ Run a simulation

        Parameters
        ----------
        tel : lsst.sims.speedObservatory.Telescope
            The Telescope instance to run the simulation for

        Notes
        -----
        Parameters controlling the operation of the Simulator are at
        the top of this file. Parameters controlling the survey
        are in config.py
        """
        self.tel = tel

        # read in downtime nights TODO use SOCS/speedObservatory for this
        downtimeNights = self._parseDowntime()

        # initialize the scheduler
        self.sched = Scheduler(telescope=self.tel, context=self)

        # make a TimeHandler to give to the CloudModel
        # TODO use sims_speedObservatory instead
        dateFormat = "%Y-%m-%d"
        startDatetime = datetime.utcfromtimestamp(config.surveyStartTime)
        timeHandler = TimeHandler(datetime.strftime(startDatetime, dateFormat))
        self.cloudModel = CloudModel(timeHandler)
        # load the cloud database
        self.cloudModel.initialize()

        # create the skyMap and display if necessary
        if trackMap:
            self.skyMap = SkyMap(telescope=self.tel, resScale=resScale)
        if showDisp:
            self.display = GraphicalMonitor(skyMap=self.skyMap, mode=graphicsMode)

        # track slew times for simplest summary statistics
        self.slewTimes = []

        # write the header to the output file if writeCsv flag is set
        if writeCsv:
            self.outFile = open("results/" + runName + "/" + runName + ".csv", "w")
            self.outFile.write("mjd,prop,ra,dec,filter\n")

        # these variables keep track of wasted time
        self.fieldsRisingWasteTime = 0
        self.earlyNightEndWasteTime = 0

        # run the survey!
        for nightNum in range(config.surveyNumNights):
            print("Night:", nightNum, end="\r")
            sys.stdout.flush()

            # skip downtime nights
            if nightNum in downtimeNights:
                continue

            # simulate the night
            self._simulateNight(nightNum)
            self.sched.notifyNightEnd()

        print("time wasted waiting for fields to rise:", self.fieldsRisingWasteTime)
        print("time wasted when sched ran out of visits:", self.earlyNightEndWasteTime)
        if writeCsv:
            self.outFile.close()

        # we're done with the simulation now
        # show summary plots if requested
        if showSummaryPlots:
            self._outputSummaryStats()

    def _simulateNight(self, nightNum):
        """ Simulate a single night

        Parameters
        ----------
        nightNum : int
            The index of the night to be scheduler. The Simulator is agnostic
            about the order in which nights are simulated, but the Scheduler
            will probably get confused if nights are not simulated sequentially
        """
        if doXTwil:
            twilStart = sky.xTwilStart(config.surveyStartTime, nightNum)
            twilEnd   = sky.xTwilEnd(config.surveyStartTime, nightNum)
        else:
            twilStart = sky.twilStart(config.surveyStartTime, nightNum)
            twilEnd   = sky.twilEnd(config.surveyStartTime, nightNum)

        # start out the simulation at the beginning of the night
        self.curTime = twilStart

        # prevI is the last value of i when we updated the display
        prevI = 0

        # prevFilter is the filter of the last visit
        (prevFilter, prevAlt, prevAz) = ('', -1, -1)

        # flag to keep track of whether the dome was just closed due to weather
        wasDomeJustClosed = False

        # loop through the visits given by the scheduler
        direction = NORTH if nightNum % 2 == 0 else SOUTH
        xTwilSched = XTwilScheduler(self, self.tel, nightNum, direction)
        #for i, visit in enumerate(self.sched.scheduleNight(nightNum)):
        for i, visit in enumerate(xTwilSched.schedule()):
            # figure out whether/how we should update the display
            perNight, deltaI = self.getUpdateRate(nightNum, i)

            # stop if the night is over
            if self.curTime >= twilEnd:
                break

            # skip forward in time if there are clouds
            deltaT = self.curTime - config.surveyStartTime
            cloudCover = self.cloudModel.get_cloud(deltaT)
            timeBeforeDomeClose = self.curTime
            while cloudCover > config.maxCloudCover and self.curTime <= twilEnd:
                self.curTime += 600
                deltaT = self.curTime - config.surveyStartTime
                cloudCover = self.cloudModel.get_cloud(deltaT)

            # if the dome closed due to clouds, get a fresh visit from the
            # scheduler by continuing
            if self.curTime > timeBeforeDomeClose:
                # let sched know that the dome closed for a while
                if not doXTwil:
                    self.sched.notifyDomeClosed(self.curTime - timeBeforeDomeClose)
                wasDomeJustClosed = True
                # continue to get a new visit that isn't stale
                continue

            if visit is not None:
                alt, az = sky.radec2altaz(visit.ra, visit.dec, self.curTime)
                # make sure this az is a valid place to look
                if alt < self.tel.minAlt or alt > self.tel.maxAlt:
                    #print("invalid alt (", np.degrees(alt), "deg) night", nightNum)
                    continue

                # figure out how far we have to slew
                slewTime = 0
                if i > 0 and not wasDomeJustClosed:
                    slewTime = self.tel.calcSlewTime(prevAlt, prevAz, prevFilter,
                                                     alt, az, visit.filter,
                                                     laxDome = config.laxDome)
                    self.curTime += slewTime
                    self.slewTimes.append(slewTime)

                # don't observe if the night is about to end
                if(self.curTime + visit.expTime +
                        config.visitOverheadTime > twilEnd):
                    break

                # notify the skyMap of the visit
                # (the time of the visit is the time after the slew is over)
                if trackMap:
                    self.skyMap.addVisit(visit, slewTime, self.curTime)

                # write the observation to the csv output if necessary
                if writeCsv:
                    assert(visit.prop == PROP_WFD or visit.prop == PROP_DD)
                    prop = "wfd" if visit.prop == PROP_WFD else "dd"
                    self.outFile.write(str(unix2mjd(self.curTime)) + "," +
                                       prop + "," +
                                       str(visit.ra) + "," +
                                       str(visit.dec) + "," +
                                       visit.filter + "\n")

                # add the exposure time of this visit to the current time
                assert(visit.expTime >= 0)
                self.curTime += visit.expTime
                self.curTime += config.visitOverheadTime

                # let the scheduler know we "carried out" this visit
                self.sched.notifyVisitComplete(visit, self.curTime)

                (prevAlt, prevAz, prevFilter) = (alt, az, visit.filter)
                wasDomeJustClosed = False
            else:
                # the visit was None, which means the scheduler has nowhere to
                # point at the moment
                self.curTime += 120
                self.fieldsRisingWasteTime += 120
                wasDomeJustClosed = True


            # now that we've added the visit (if there was one),
            # update the display
            if showDisp and not perNight and i - prevI >= deltaI:
                self.display.updateDisplay(self.skyMap, self.curTime)
                prevI = i
                # save each frame if the saveMovie flag is set
                if saveMovie:
                    self.display.saveFrame(movieDir + "/%04d%04d.png" % (nightNum, i))

        # the night is over

        # calculate how much time we waste at the end of the night
        # curTime could be > twilEnd if we added 600 sec for clouds
        self.earlyNightEndWasteTime += max(twilEnd - self.curTime, 0)
        # show the sky scroll by during the wasted time
        if showDisp and not perNight:
            while self.curTime < twilEnd:
                self.display.updateDisplay(self.skyMap, self.curTime)
                self.curTime += 30

        # update the display with the completed night if perNight is true
        if showDisp and perNight:
            self.display.updateDisplay(self.skyMap, self.curTime)
            if saveMovie:
                self.display.saveFrame(movieDir + "/%07d.png" % i)

        # clear the skyMap so the next updateDisplay won't have tonight's visits
        if showDisp and clearDisplayNightly:
            self.skyMap.clear()
        if 0 < len(self.slewTimes) < 1:
            print("S" if direction == SOUTH else "N", end="")
            print(nightNum, end=":")
            print("#/mean slew:", len(self.slewTimes), np.mean(self.slewTimes))
        self.slewTimes = []

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

        print("avg slew time", np.mean(self.slewTimes), "seconds")
        print("median slew time", np.median(self.slewTimes), "seconds")
        plotter.slewHist()
        
        sortedTimes = np.sort(self.slewTimes)
        cum = np.cumsum(sortedTimes)
        print("total cumulative slew time: ", cum[-1])
        print("rank @ half total cum / # slews", np.searchsorted(cum, cum[-1]/2) / len(cum))
        plotter.slewRankCum()
        
        plotter.revisitHist()
        
        plotter.dAirmassCum()
        plotter.dAirmassContour()
        plotter.zenithAngleContour()
        plotter.airmassHist()

        plotter.show()

    @staticmethod
    def getUpdateRate(nightNum, i):
        """ Controls how often the display is updated

        Parameters
        ----------
        nightNum : int
            The index of the night being simulated
        i : int
            The index of the visit within the night

        Returns
        -------
        perNight : bool
            Whether or not to update the display nightly
        deltaI : int
            If perNight is False, how often to update the display

        Notes
        -----
        This method is only useful if you want a variable update rate --
        i.e. if you want to make a movie of a simulation that changes speed
        over time. What's currently implemented will make the video speed
        up from one frame per visit up to one frame per night.
        """
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
    def _parseDowntime(self, schedFileName="schedDown.conf",
                            unschedFileName="unschedDown.conf"):
        """ Gets the downtime nights from the configuration files

        Parameters
        ----------
        schedFileName : string
            The name of the file containing the scheduled downtime
        unschedFileName : string
            The name of the file containing the unscheduled downtime

        Returns
        -------
        A set containing all night indices that are down (whether scheduled
        or unscheduled)
        """
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
    """ Runs a simulation with the default parameters """
    sim = Simulator()
    # use a telescope with default parameters
    tel = Telescope()

    # make the results directory if necessary
    if not os.path.isdir("results/" + runName):
        os.mkdir("results/" + runName)
    return sim.run(tel)

if __name__ == "__main__":
    runDefaultSim()
