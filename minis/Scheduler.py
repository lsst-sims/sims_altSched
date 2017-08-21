from __future__ import division

import numpy as np
import Config
from Config import NORTH, SOUTH, EAST, WEST, SOUTHEAST
from Visit import PROP_WFD, PROP_DD
from lsst.sims.speedObservatory import sky
from lsst.sims.speedObservatory import Telescope
from Visit import Visit
from Visit import VisitPair
import Utils
from NightScheduler import NightScheduler

from matplotlib import pyplot as plt
import copy
import csv

class Scheduler:
    def __init__(self, telescope, context):
        self.telescope = telescope
        self.context = context
        # or maybe load from .pkl
        self.makeupVPs = set()

        # keep track of how long we spend slewing each night
        # so we can update our estimate of the average slew time each night
        self.tonightsSlewTimes = []



        # keep track of how many visits have been executed in each direction
        self.SVisitsComplete = 0
        self.NVisitsComplete = 0
        self.EVisitsComplete = 0

    def scheduleNight(self, nightNum):
        # decide which way to point tonight
        NCoverage = self.NVisitsComplete / Utils.areaInDir(NORTH)
        SCoverage = self.SVisitsComplete / Utils.areaInDir(SOUTH)
        ECoverage = self.EVisitsComplete / Utils.areaInDir(EAST)
        SECoverage = ((self.SVisitsComplete + self.EVisitsComplete) /
                      Utils.areaInDir(SOUTHEAST))

        if NCoverage < SECoverage:
            self.nightDirection = NORTH
        else:
            self.nightDirection = SOUTHEAST

        # reset the slew times array
        self.tonightsSlewTimes = []
        prevAlt = prevAz = None
        prevFilter = self.telescope.filters[0]

        # return each visit prescribed by tonight's NightScheduler
        self.nightScheduler = NightScheduler(self.telescope, nightNum,
                                             self.nightDirection, self.makeupVPs)
        prevTime = None
        for visit in self.nightScheduler.schedule():
            time = self.context.time()
            alt, az = sky.radec2altaz(visit.ra, visit.dec, self.context.time())
            if alt < self.telescope.minAlt:
                # East is +pi/2, so if the field has az < pi, it is rising
                # and if az > pi then setting
                if az >= np.pi:
                    # this field is setting, so skip it
                    continue
                else:
                    # this field is rising, so wait a while until it's
                    # visible
                    while alt < self.telescope.minAlt:
                        # if we yield None the simulator (or the world) will
                        # progress time for us
                        yield None
                        alt, az = sky.radec2altaz(visit.ra, visit.dec,
                                                  self.context.time())
                        prevAlt = prevAz = None
            if prevAlt is not None:
                # Don't change laxDome param without changing in Simulator too
                slewTime = self.telescope.calcSlewTime(prevAlt, prevAz, prevFilter,
                                                       alt, az, visit.filter,
                                                       laxDome = Config.laxDome)
                self.tonightsSlewTimes.append(slewTime)
            prevAlt = alt
            prevAz = az
            prevFilter = visit.filter
            prevTime = time
            yield visit

    def notifyVisitComplete(self, visit, time):
        if not isinstance(visit, Visit):
            raise TypeError("must pass in Visit() instance")

        if visit.isComplete:
            raise RuntimeError("visit was completed twice")

        visit.isComplete = True
        visit.timeOfCompletion = time

        # check if this visit pair is complete
        visitPair = visit.visitPair
        if visitPair == None:
            # this visit is not part of a visit pair, so it must be
            # a DD visit
            assert(visit.prop == PROP_DD)
            # we don't currently keep track of which DD visits are
            # actually carried out, so just return
            return

        # visits with a non-None visit pair are WFD visits
        assert(visit.prop == PROP_WFD)
        if visitPair.visit1.isComplete and visitPair.visit2.isComplete:
            # remove the visit pair from either tonights mini survey
            # or from the makeup VPs (whichever it happens to be in)
            # (set.discard removes if the element exists)
            self.makeupVPs.discard(visitPair)
            self.nightScheduler.notifyVisitPairComplete(visitPair)

            # keep track of our coverage in each direction
            if Utils.directionOfDec(visitPair.dec) == NORTH:
                self.NVisitsComplete += 1
            elif Utils.directionOfDec(visitPair.dec) == SOUTH:
                self.SVisitsComplete += 1
            elif Utils.directionOfDec(visitPair.dec) == EAST:
                self.EVisitsComplete += 1
            else:
                raise RuntimeError("Completed visit " + str(visit) + \
                                   " is in unknown direction")

    def notifyNightEnd(self):
        if len(self.tonightsSlewTimes) == 0:
            return
        avgSlew = np.mean(self.tonightsSlewTimes)
        newMakeups = self.nightScheduler.notifyNightEnd(avgSlew)
        self.makeupVPs.update(newMakeups)

    def notifyDomeClosed(self, timeClosed):
        self.nightScheduler.notifyDomeClosed(timeClosed)
