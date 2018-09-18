from __future__ import division, print_function

from config import NORTH, SOUTH, EAST, WEST, SOUTHEAST
import config
from minis.MiniSurvey import MiniSurvey
from Visit import Visit
#from Visit import PROP_XTWIL
from minis import filtersequence
import numpy as np

from lsst.sims.speedObservatory import sky

class XTwilScheduler:
    """ Class that schedules the extended twilight time of a night

    """
    def __init__(self, context, telescope, nightNum, direction):
        self.context = context
        self.telescope = telescope
        self.nightNum = nightNum
        self.direction = direction

        # width of a scan in units of RA
        if direction == NORTH:
            #scanWidth = 2.8 * self.telescope.fovWidth # default zenithBuffer
            scanWidth = 2 * self.telescope.fovWidth
            self.numScans = 8

            self.minDec = np.radians(-15)
            self.maxDec = np.radians(5)

            self.minDec = np.radians(5)
            self.maxDec = np.radians(22)
        elif direction == SOUTH:
            #scanWidth = 4.5 * self.telescope.fovWidth # default zenithBuffer
            scanWidth = 3.5 * self.telescope.fovWidth
            self.numScans = 4

            self.minDec = np.radians(-66)
            self.maxDec = np.radians(-40)

            self.minDec = np.radians(-90)
            self.maxDec = np.radians(-66)

        # we have time to do about 2 scans in each extended
        # twilight period (sunset and sunrise)
        #self.numScans = 4

        sst = config.surveyStartTime
        self.eveningStartTime = sky.xTwilStart(sst, self.nightNum)
        self.eveningEndTime   = sky.twilStart( sst, self.nightNum)
        self.morningStartTime = sky.twilEnd(   sst, self.nightNum)
        self.morningEndTime   = sky.xTwilEnd(  sst, self.nightNum)

        self.eveningStartRa = sky.raOfMeridian(self.eveningStartTime) + self.telescope.fovWidth * 4
        self.eveningEndRa   = self.eveningStartRa + self.numScans * scanWidth
        self.morningEndRa   = sky.raOfMeridian(self.morningEndTime) - self.telescope.fovWidth * 4
        self.morningStartRa = self.morningEndRa - self.numScans * scanWidth

        self.eveningStartRa %= (2*np.pi)
        self.eveningEndRa   %= (2*np.pi)
        self.morningStartRa %= (2*np.pi)
        self.morningEndRa   %= (2*np.pi)

        self.eveningVPs = MiniSurvey.newMiniSurvey(
                    self.telescope, self.eveningStartRa, self.eveningEndRa, direction,
                    minDec=self.minDec, maxDec=self.maxDec
                )
        self.morningVPs = MiniSurvey.newMiniSurvey(
                    self.telescope, self.morningStartRa, self.morningEndRa, direction,
                    minDec=self.minDec, maxDec=self.maxDec
                )

    def schedule(self):
        """ Schedule the extended twilight time

        """
        if self.direction == SOUTH:
            for visit in self._scheduleSPCap():
                yield visit
            return

        morningVPs = np.array(list(self.morningVPs))
        eveningVPs = np.array(list(self.eveningVPs))

        i = 0
        for visit in self._scheduleVPs(eveningVPs, self.eveningStartTime, self.eveningEndTime):
            i += 1
            yield visit
        if i == eveningVPs.size:
            #print("ran out of evening visits in North on night {}".format(self.nightNum))
            pass
        for visit in self._scheduleVPs(morningVPs, self.morningStartTime, self.morningEndTime):
            yield visit

        #print("ran out of morning visits in North on night {}".format(self.nightNum))
        return

        for VPs, label in zip([self.eveningVPs, self.morningVPs],
                              ["evening", "morning"]):
            VPs = np.array(list(VPs))

            # divide into two scans, then execute one up and one down
            scans = self._getScans(VPs)
            scanDirs = self._getScanDirs()


            for scan, scanDir in zip(scans, scanDirs):
                if scanDir == NORTH:
                    sortedScan = sorted(scan, key=lambda v: v.dec)
                elif scanDir == SOUTH:
                    sortedScan = sorted(scan, key=lambda v: -1*v.dec)
                
                # we only schedule one visit from each visit pair --
                # no revisits for the extended twilight time
                for VP in sortedScan:
                    eveningOver    = self.context.time() > self.eveningEndTime
                    morningEnded   = self.context.time() > self.morningEndTime
                    if label == "evening" and eveningOver:
                        break
                    while label == "morning" and (self.context.time() <
                                                  self.morningStartTime):
                        yield None
                    if label == "morning" and morningEnded:
                        break
                    VP.visit1.expTime = config.xTwilExpTime
                    VP.visit1.filter = "g"
                    yield VP.visit1
            if ((label == "evening" and self.context.time() < self.eveningEndTime) or 
                (label == "morning" and self.context.time() < self.morningEndTime)):
                print("{}: ran out of visits in {}".format(self.nightNum, label))

    def _scheduleVPs(self, VPs, startTime, endTime):
        scans = self._getScans(VPs)
        scanDirs = self._getScanDirs()

        for scan, scanDir in zip(scans, scanDirs):
            if scanDir == NORTH:
                sortedScan = sorted(scan, key=lambda v: v.dec)
            elif scanDir == SOUTH:
                sortedScan = sorted(scan, key=lambda v: -1*v.dec)
            
            # we only schedule one visit from each visit pair --
            # no revisits for the extended twilight time
            for VP in sortedScan:
                while self.context.time() < startTime:
                    # not ready to schedule yet, yield None until
                    # enough time passes
                    yield None
                if self.context.time() > endTime:
                    # out of time, so return
                    return

                VP.visit1.expTime = config.xTwilExpTime
                VP.visit1.filter = "g"
                yield VP.visit1

    def _scheduleSPCap(self):
        # get VPs around the south pole
        VPs = MiniSurvey.newMiniSurvey(
                self.telescope,
                #(self.eveningStartRa - np.radians(120)) % (2*np.pi),
                #(self.morningEndRa   + np.radians(120)) % (2*np.pi),
                0, 2*np.pi,
                SOUTH, minDec=self.minDec, maxDec=self.maxDec
            )
        VPs = np.array(list(VPs))
        try:
            assert(VPs.size > 0)
        except AssertionError as e:
            print("ra range: {}-{}".format((self.eveningStartRa - np.radians(120)) % (2*np.pi), (self.morningEndRa   + np.radians(120)) % (2*np.pi)))
            raise e

        # supriously rotate the VPs so that we can easily turn into scans
        # (normally we divide into scans by RA but this doesn't work at the poles)
        for VP in VPs:
            newRa, newDec = self._rotate(VP.ra, VP.dec, np.pi/2)
            # note this leaves the ra/dec of the VP's visits (on purpose), which is
            # terrible coding. Sorry
            VP.ra = newRa
            VP.dec = newDec

        i = 0
        for visit in self._scheduleVPs(VPs, self.eveningStartTime, self.eveningEndTime):
            assert(visit is not None or i > 0)
            if visit is not None:
                i += 1
            yield visit
        if i == VPs.size:
            print("finished whole night in evening")
            return
        remainingVPs = np.array([VP for VP in VPs if not VP.visit1.isComplete])
        for visit in self._scheduleVPs(remainingVPs, self.morningStartTime, self.morningEndTime):
            yield visit
        #print("ran out of visits in the South on night {}".format(self.nightNum))

    def _getScanDirs(self):
        firstDirection = self.direction
        scanDirs = [firstDirection if i % 2 == 0
                    else NORTH if firstDirection == SOUTH else SOUTH
                    for i in range(self.numScans)]
        return scanDirs


    def _getScans(self, VPs):
        ras = np.array([VP.ra for VP in VPs])
        # ra range will always be less than 180 degrees, so this
        # means ras wraps around 2pi
        if ras.max() - ras.min() > np.pi:
            ras = (ras - np.pi) % (2*np.pi)

        scanBoundaryRas = [np.percentile(ras, i * 100 / self.numScans) 
                           for i in range(1, self.numScans)]
        scanBoundaryRas = [ras.min() - 0.1] + scanBoundaryRas + [ras.max() + 0.1]
        scans = [VPs[(ras >= scanBoundaryRas[i]) & (ras < scanBoundaryRas[i+1])]
                 for i in range(self.numScans)]
        return scans



    def _rotate(self, ra, dec, rotAmt):
        # rotate 90 degrees about the x axis
        sin = np.sin
        cos = np.cos
        phi = ra
        theta = dec
        newRa = np.arctan2(-sin(rotAmt)*sin(theta) + cos(rotAmt)*cos(theta)*sin(phi),
                            cos(theta)*cos(phi))
        newRa %= 2 * np.pi
        newDec = np.arcsin(cos(rotAmt)*sin(theta) + sin(rotAmt)*cos(theta)*sin(phi))
        return (newRa, newDec)

    def notifyVisitComplete(self, visit, time):
        self.think("lol")
        return

    def notifyDomeClosed(self, closeTime):
        self.think("hah")
        return
    
    def think(self, msg):
        if msg == "segfault":
            raise Exception()
        return
