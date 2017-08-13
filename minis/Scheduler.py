from __future__ import division
import numpy as np
import Config
from Config import NORTH, SOUTH, EAST, WEST, SOUTHEAST
from Visit import PROP_WFD, PROP_DD
from lsst.sims.speedObservatory import sky
from lsst.sims.speedObservatory import Telescope
from MiniSurvey import MiniSurvey
import filtersequence
from Visit import Visit
from Visit import VisitPair
import Utils

from matplotlib import pyplot as plt
import copy
import csv


# add a new mini-survey whenever we dip below this much times
# the expected number of visits needed in a night
NUM_VISITS_PADDING = 1

class Scheduler:
    def __init__(self, telescope, context):
        self.telescope = telescope
        self.context = context
        # or maybe load from .pkl
        self.makeupVPs = set()

        # these estimates get updated at the end of each night
        # (currently estAvgExpTime is not updated)
        self.NEstAvgSlewTime = 10 # seconds
        self.SEEstAvgSlewTime = 10 # seconds
        self.estAvgExpTime = 30 #seconds

        # keep track of how long we spend slewing each night
        # so we can update our estimate of the average slew time each night
        self.curNightSlewTimes = []

        # int_{dec1}^{dec2} 2\pi r dr, where r=\cos\theta
        # = 2\pi (\sin(dec2) - \sin(dec1)) 
        buf = Config.zenithBuffer
        self.areaInNorth = 2*np.pi*(np.sin(Config.maxDec) - \
                                    np.sin(self.telescope.latitude + buf))

        # similar for south and east
        self.areaInSouth = 2*np.pi*(np.sin(self.telescope.latitude - buf) -
                                    np.sin(Config.minDec))
        self.areaInEast  = 2*np.pi*(np.sin(self.telescope.latitude + buf) -
                                    np.sin(self.telescope.latitude - buf))

        # keep track of how many visits have been executed in each direction
        self.SVisitsComplete = 0
        self.NVisitsComplete = 0
        self.EVisitsComplete = 0

        # read in DD field centers
        self.ddFieldCenters = []
        with open("ddFields.csv", "r") as ddFile:
            fields = csv.DictReader(ddFile)
            for field in fields:
                ra = float(field["ra"])
                dec = float(field["dec"])
                self.ddFieldCenters.append((ra, dec))


    def _getFilters(self, nightNum, direction, nScans):
        return filtersequence.maxFChanges(nightNum, direction, nScans)

    def scheduleNight(self, nightNum):
        # decide which way to point tonight
        NCoverage = self.NVisitsComplete / self.areaInNorth
        SCoverage = self.SVisitsComplete / self.areaInSouth
        ECoverage = self.EVisitsComplete / self.areaInEast
        SECoverage = ((self.SVisitsComplete + self.EVisitsComplete) /
                      (self.areaInSouth     + self.areaInEast))
        #print "NCoverage", NCoverage
        #print "SCoverage", SCoverage
        #print "ECoverage", ECoverage
        #print "SECoverage", SECoverage
        if NCoverage < SECoverage:
            self.nightDirection = NORTH
        else:
            self.nightDirection = SOUTHEAST

        # reset the slew times array
        self.curNightSlewTimes = []
        prevAlt = prevAz = None
        prevFilter = self.telescope.filters[0]

        # return each visit prescribed by _scheduleNight()
        prevTime = None
        for visit in self._scheduleNight(nightNum):
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
                                                       laxDome = True)
                self.curNightSlewTimes.append(slewTime)
            prevAlt = alt
            prevAz = az
            prevFilter = visit.filter
            prevTime = time
            yield visit

    # TODO confusing why we have scheduleNight() and _scheduleNight()
    def _scheduleNight(self, nightNum):
        # this will hold the VPs in tonight's mini survey
        self.tonightsMini = None

        # figure out what DD visit we're going to do tonight
        self.DDVisit = None
        nightStart = sky.nightStart(Config.surveyStartTime, nightNum)
        nightEnd = sky.nightEnd(Config.surveyStartTime, nightNum)
        minRa = sky.raOfMeridian(nightStart)
        maxRa = sky.raOfMeridian(nightEnd)
        for ra, dec in self.ddFieldCenters:
            # check if the night will end before we could finish
            # the DD exposure
            # raBuffer is how far the sky will move during the exposure
            raBuffer = Config.DDExpTime / (24*3600) * 2*np.pi

            DDDirection = self._directionOfDec(dec)

            if DDDirection == NORTH:
                inRaRange = Utils.isRaInRange(ra, (minRa, maxRa - raBuffer))
            elif DDDirection == SOUTH:
                inRaRange = Utils.isRaInRange(ra, (minRa, maxRa - raBuffer))
                DDDirection = SOUTHEAST
            elif DDDirection == EAST:
                # DD visits in the East need to be in the offset
                # allowable RA range
                raRange = (minRa + Config.zenithBufferOffset + raBuffer,
                           maxRa + Config.zenithBufferOffset - raBuffer)
                inRaRange = Utils.isRaInRange(ra, raRange)
                DDDirection = SOUTHEAST

            # the DD ra must be less than maxRa - raBuffer so that
            # the night doesn't end while the DD exposure is happening
            if inRaRange and DDDirection == self.nightDirection:
                self.DDVisit = Visit(PROP_DD, None, ra, dec, 0, Config.DDExpTime, 'u')
                break
        # XXX don't return DD
        #self.DDVisit = None

        # helper for later to check if a DD visit falls within the min/max
        # RA range of a revisit group
        def isDDVisitInRevisitGroup(revisitGroup):
            scan1 = revisitGroup["scan1"]
            scan2 = revisitGroup["scan2"]
            if len(scan1) == 0 and len(scan2) == 0:
                    return False
            ras = np.array([vp.ra for vp in revisitGroup["scan1"]] +
                           [vp.ra for vp in revisitGroup["scan2"]])
            if ras.max() - ras.min() > np.pi:
                # scans are never > 180 degrees in ra, so if true, the scan
                # must cross the 2pi boundary
                maxRa = ((ras + np.pi) % (2*np.pi)).max() - np.pi
                minRa = ((ras + np.pi) % (2*np.pi)).min() - np.pi
            else:
                maxRa = ras.max()
                minRa = ras.min()
            return Utils.isRaInRange(self.DDVisit.ra, (minRa, maxRa))


        # figure out how many visits we'll probably do tonight
        nightLength = sky.nightLength(Config.surveyStartTime, nightNum)
        # the period of exposure taking is the sum of the exposure time,
        # the visit overhead time (shutter/intermediate readout),
        # and the slew time
        if self.nightDirection == NORTH:
            t = self.NEstAvgSlewTime
        elif self.nightDirection == SOUTHEAST:
            t = self.SEEstAvgSlewTime
        t += self.estAvgExpTime + Config.visitOverheadTime
        if self.DDVisit is None:
            # we won't do a DD visit tonight
            expectedNumVisits = int(nightLength / t)
        else:
            # we do plan to do a DD visit tonight, so the amount of time
            # available for WFD exposures is nightLength - DDexpTime
            # TODO this doesn't include slew time to/from the DD field
            expectedNumVisits = int((nightLength - self.DDVisit.expTime) / t)

        # TODO it might be greater than numVisits / 2 if we expect to
        # lose revisits due to weather. May model this at some point
        # (VP = Visit Pair)
        expectedNumVPs = int(expectedNumVisits / 2)

        """ here is the scanning pattern:

        the sky moves to the right over time so that the
        vertical scans are all right around the meridian,
        which in this diagram is vertical


           scan2       scan1        scan0
        -------------------------------------
        |           |       ||<=|=========  |
        |           |  ||<==||==|====   /\  |
        |           |  ||   ||  |  /\   ||  |
        |           |  ||   ||  |  ||   ||  |
        |           |  ||   ||  |  ||   ||  |
        |           |  ||   ||  |  ||   ||  |
        |           |  ||   ||  |  ||   ||  |
        |      etc. |  ||   ||  |  ||   ||  |
        |       /\  |  ||   ||  |  ||   ||  |
        |       ||  |  ||   ||  |  ||   ||  |
        |       ||  |  \/   \/  |  ||   ||  |
        |       ||<=|====   ====|=>||  start|
        -------------------------------------

        """

        if self.nightDirection == NORTH:
            makeupVPs = self._getMakeupVPs(nightNum, NORTH)
            if len(makeupVPs) >= expectedNumVPs * NUM_VISITS_PADDING:
                # tonight will be a purely makeup night
                self.tonightsMini = set()
                tonightsVPs = np.random.choice(list(makeupVPs),
                                               size=int(expectedNumVPs),
                                               replace=False)
            else:
                # we start with one new minisurvey for tonight and then
                # fill in the remaining needed visits with makeup VPs
                self.tonightsMini = self._getNewMiniSurvey(nightNum, NORTH)

                # add to makeupVPs if we still don't have enough visit pairs
                while(len(self.tonightsMini) + len(makeupVPs) <
                        expectedNumVPs * NUM_VISITS_PADDING):
                    self.makeupVPs.update(self._getNewMiniSurvey(nightNum, NORTH))
                    makeupVPs = self._getMakeupVPs(nightNum, NORTH)
                tonightsVPs = list(self.tonightsMini)
                if len(self.tonightsMini) < int(expectedNumVPs):
                    tonightsVPs.extend(
                        np.random.choice(list(makeupVPs),
                                         size=int(expectedNumVPs) - len(self.tonightsMini),
                                         replace=False)
                    )

            # choose which visitPairs we're going to execute tonight
            # TODO could prioritize doing stale visitPairs
            scans = self._generateScans(nightNum, tonightsVPs, NORTH)
            filters = self._getFilters(nightNum, NORTH, len(scans)*2)

            # a revisitGroup is a pair of scans with info on how to
            # execute them
            revisitGroups = []
            for i in range(0, len(scans), 2):
                revisitGroup = {"scan1": scans.pop(0),
                                "scan2": scans.pop(0),
                                "scanDir1": NORTH,
                                "scanDir2": SOUTH,
                                "filter1": filters.pop(0),
                                "filter2": filters.pop(0),
                                "filter3": filters.pop(0),
                                "filter4": filters.pop(0)}
                revisitGroups.append(revisitGroup)

            # we should have used up all the scans and filters
            assert(len(filters) == 0)
            assert(len(scans) == 0)

        elif self.nightDirection == SOUTHEAST:
            areaInSouthEast = self.areaInSouth + self.areaInEast
            SExpectedNumVPs = expectedNumVPs * self.areaInSouth / areaInSouthEast
            EExpectedNumVPs = expectedNumVPs * self.areaInEast / areaInSouthEast

            SMakeupVPs = self._getMakeupVPs(nightNum, SOUTH)
            EMakeupVPs = self._getMakeupVPs(nightNum, EAST)
            if len(SMakeupVPs) >= SExpectedNumVPs * NUM_VISITS_PADDING and \
               len(EMakeupVPs) >= EExpectedNumVPs * NUM_VISITS_PADDING:
                # tonight is a makeup  only night
                self.tonightsMini = set()

                STonightsVPs = np.random.choice(list(SMakeupVPs),
                                                size=int(SExpectedNumVPs),
                                                replace=False)
                ETonightsVPs = np.random.choice(list(EMakeupVPs),
                                                size=int(EExpectedNumVPs),
                                                replace=False)
            else:
                # make a new mini for the south and east, then fill in each
                # direction with makeup VPs
                SMini = self._getNewMiniSurvey(nightNum, SOUTH)
                EMini = self._getNewMiniSurvey(nightNum, EAST)
                self.tonightsMini = SMini | EMini

                # if needed, add minis to S and E makeupVPs
                while(len(SMini) + len(SMakeupVPs) <
                        SExpectedNumVPs * NUM_VISITS_PADDING):
                    self.makeupVPs.update(self._getNewMiniSurvey(nightNum, SOUTH))
                    SMakeupVPs = self._getMakeupVPs(nightNum, SOUTH)

                while(len(EMini) + len(EMakeupVPs) <
                        EExpectedNumVPs * NUM_VISITS_PADDING):
                    self.makeupVPs.update(self._getNewMiniSurvey(nightNum, EAST))
                    EMakeupVPs = self._getMakeupVPs(nightNum, EAST)
 
                # tonights VPs in the South/East consist of tonights mini in
                # that direction plus a random selection of makeups from that
                # direction
                STonightsVPs = list(SMini)
                ETonightsVPs = list(EMini)
                if len(SMini) < int(SExpectedNumVPs):
                    STonightsVPs.extend(
                        np.random.choice(list(SMakeupVPs),
                                         size=int(SExpectedNumVPs) - len(SMini),
                                         replace=False)
                    )
                if len(EMini) < int(EExpectedNumVPs):
                    ETonightsVPs.extend(
                        np.random.choice(list(EMakeupVPs),
                                         size=int(EExpectedNumVPs) - len(EMini),
                                         replace=False)
                    )
            #print "S tonight", len(STonightsVPs) #/ self.areaInSouth
            #print "E tonight", len(ETonightsVPs) #/ self.areaInEast

            SScans = self._generateScans(nightNum, STonightsVPs, SOUTH)
            EScans = self._generateScans(nightNum, ETonightsVPs, EAST)

            # TODO confusing use of the word scan to mean both a single vertical
            # strip of visits (that gets executed twice) and a single execution
            # of a verticle strip.
            nSScans = len(SScans)
            nScans = (len(SScans) + len(EScans)) * 2
            filters = self._getFilters(nightNum, SOUTHEAST, nScans)

            #print "S scans", len([v for s in SScans for v in s]) / self.areaInSouth
            #print "E scans", len([v for s in EScans for v in s]) / self.areaInEast

            # do as many southern pairs as we can before we have to
            # catch the E fields before they hit the zenith avoidance zone
            raOfZenith, _ = sky.altaz2radec(np.pi/2, 0,
                                            self.context.time())
            cutoffRa = raOfZenith + Config.zenithBuffer

            EScansRas = [np.array([v.ra for v in scan]) for scan in EScans]
            EMinRas = np.zeros(len(EScansRas))
            for i, ras in enumerate(EScansRas):
                if len(ras) == 0:
                    # TODO
                    #print "len(ras)=0 :(. nightNum =", nightNum
                    #print "nightLen", sky.nightLength(Config.surveyStartTime, nightNum) / 3600
                    #print "EScans len", [len(scan) for scan in EScans]
                    #print "SScans len", [len(scan) for scan in SScans]
                    #print "avgVisitTime", self.estAvgExpTime + self.SEEstAvgSlewTime
                    #raise RuntimeError("ras len 0, i=%d" % i)

                    # this will cause this scan to be scheduled ASAP
                    # since the scheduler will think there's no time left
                    # before it hits zenith
                    # But this is fine since there are no visits in this scan
                    # anyway.
                    EMinRas[i] = cutoffRa
                    continue
                if ras.max() - ras.min() > np.pi:
                    minRa = ((ras + np.pi) % (2*np.pi)).min() - np.pi
                else:
                    minRa = ras.min()
                EMinRas[i] = minRa
                #EMinRas = np.array([min(ras) for ras in EScansRas])
            numECols = int(len(EScans) / 2)
            # EMinGroupRa is an array of numECols numbers: the min RA from
            # scans 1/2, from scans 3/4, 5/6, etc.
            EMinGroupRa = EMinRas.reshape(numECols,2).T.min(axis=0)

            timesLeft = ((EMinGroupRa - cutoffRa) % (2*np.pi)) * (3600*24) / (2*np.pi)

            avgVisitTime = (self.estAvgExpTime + self.SEEstAvgSlewTime +
                            Config.visitOverheadTime)
            EScanTimes = np.array(map(len, EScans)) * avgVisitTime

            # this is the time available between when we start a pair of East
            # scans and when we have to be done with those scans. It is NOT
            # the amount of time it would take to complete each pair of East
            # scans (that would be 2 * (EScanTimes[::2] + EScanTimes[1::2]))
            # but rather the time to do the South scan once and the North scan
            # twice, since at that point we're headed back East again so are no
            # longer worried about hitting the zenith avoid zone
            execTimes = 2 * EScanTimes[::2] + EScanTimes[1::2]

            def printDebug():
                print "nightNum", nightNum
                print "numECols", numECols
                print "EMinGroupRa", EMinGroupRa
                print "raOfZenith", raOfZenith
                print "cutoffRa", cutoffRa
                print "timesLeft", timesLeft, timesLeft / 3600
                print "avgVisitTime", avgVisitTime
                print "EScanTimes", EScanTimes, EScanTimes / 3600, np.cumsum(EScanTimes/3600)
                print "execTimes", execTimes
                SScanTimes = np.array(map(len, SScans)) * avgVisitTime
                print "SScanTimes", SScanTimes, SScanTimes / 3600, np.cumsum(SScanTimes/3600)

            if nightNum > 320e10:
                printDebug()

            # now combine the South and East scans together in a good order
            # into the scans, filters, and scanDirs arrays

            # helper to make a revisitGroup -- a pair of scans with info
            # about how to execute them
            def newRevisitGroup(direction):
                assert(direction == EAST or direction == SOUTH)
                scans = EScans if direction == EAST else SScans
                revisitGroup = {"scan1": scans.pop(0),
                                "scan2": scans.pop(0),
                                "scanDir1": WEST if direction == EAST else SOUTH,
                                "scanDir2": EAST if direction == EAST else NORTH,
                                "filter1": filters.pop(0),
                                "filter2": filters.pop(0),
                                "filter3": filters.pop(0),
                                "filter4": filters.pop(0)}
                return revisitGroup


            # a revisitGroup is a pair of scans with info about how to
            # execute them
            revisitGroups = []
            cumTime = 0
            # j keeps track of which East scan we should add next
            j = 0
            # continue adding Southern scans until we need to add an East scan
            # (i.e. if we added another South scan before the next East scan,
            #  we wouldn't have time to complete the East scan before it hit
            #  the zenith avoidance zone)
            for i in range(0, nSScans, 2):
                # consider whether to add Southern scan i

                # visits get popped off the front of SScans so we always
                # consider positions 0 and 1
                numNewVisits = len(SScans[0]) + len(SScans[1])
                # this is how long it would take us to complete the next
                # Southern scan
                cumTime += 2 * numNewVisits * avgVisitTime
                paddedTime = cumTime + 0.1 * numNewVisits * avgVisitTime

                # add an East scan if we've run out of time
                # before the scan hits the zenith avoid zone
                if j < numECols and paddedTime > timesLeft[j] - execTimes[j]:
                    ERevisitGroup = newRevisitGroup(EAST)
                    revisitGroups.append(ERevisitGroup)
                    ENumNewVisits = (len(ERevisitGroup["scan1"]) +
                                     len(ERevisitGroup["scan2"]))
                    cumTime += 2 * ENumNewVisits * avgVisitTime

                    # add to cumTime if we know we'll schedule the DD visit
                    # after this revisit group
                    if (self.DDVisit is not None and
                            isDDVisitInRevisitGroup(ERevisitGroup)):
                        # TODO this doesn't include the slew time to/from
                        # the DD field
                        cumTime += self.DDVisit.expTime
                    j += 1
                # now add the South scan. Note this assumes that we never need
                # to add two East scans in a row in order to keep out of the
                # zenith avoidance region, which seems like a reasonable assumption
                # TODO but might be false in the case of extreme weather causing
                # a buildup in one southern scan
                SRevisitGroup = newRevisitGroup(SOUTH)

                # add to cumTime if we know we'll schedule the DD visit after
                # this revisit group
                if (self.DDVisit is not None and
                        isDDVisitInRevisitGroup(SRevisitGroup)):
                    # TODO this ignores the slew time to/from the DD field
                    cumTime += self.DDVisit.expTime
                revisitGroups.append(SRevisitGroup)

            # we should have added all the East scans by now excepting perhaps one
            # but sometimes if a bunch of visits pile up due to weather, we'll
            # need to add more than one E col at the end, hence the while
            # instead of if
            while j < numECols:
                revisitGroups.append(newRevisitGroup(EAST))
                j += 1

            # we should have put all the filters and scans into revisit groups
            assert(len(filters) == 0)
            assert(len(EScans) == 0)
            assert(len(SScans) == 0)

        # now schedule each revisit group
        for rg in revisitGroups:
            # perform scan1 and scan2 twice in a row but use
            # potentially-different filters the second time through
            for scan, scanDir, filt in [(rg["scan1"], rg["scanDir1"], rg["filter1"]),
                                        (rg["scan2"], rg["scanDir2"], rg["filter2"]),
                                        (rg["scan1"], rg["scanDir1"], rg["filter3"]),
                                        (rg["scan2"], rg["scanDir2"], rg["filter4"])]:
                if scanDir == NORTH:
                    sortedScan = sorted(scan, key=lambda v: v.dec)
                elif scanDir == SOUTH:
                    sortedScan = sorted(scan, key=lambda v: -1*v.dec)
                elif scanDir == EAST:
                    (raMin, raMax) = self._getNightRaRange(nightNum)
                    # subtract raMin so we never cross 2\pi (since a scan
                    # is never more than 2\pi radians long in RA)
                    sortedScan = sorted(scan, key=lambda v: (v.ra - raMin) % (2*np.pi))
                elif scanDir == WEST:
                    (raMin, raMax) = self._getNightRaRange(nightNum)
                    sortedScan = sorted(scan, key=lambda v: -1 * (v.ra - raMin) % (2*np.pi))
                else:
                    raise RuntimeError("invalid direction " + str(scanDir))

                for visit in self._schedulePath(sortedScan, nightNum, filt):
                    yield visit

            # figure out whether we need to add a DD field between
            # revisitGroups. If self.DDVisit is None, we've
            # already completed the visit
            if self.DDVisit is None:
                continue

            # now check if we're at the correct point during the night to
            # execute the DD visit

            # TODO this will miss the DD if it happens to fall exactly between
            # the maxRa of one revisitGroup and the minRa of the next
            # we should really calculate the minimum ra_{meridian} - ra_{dd}
            # over all times between revisit groups and execute it there
            if self.DDVisit is not None and isDDVisitInRevisitGroup(rg):
                yield self.DDVisit
                # we don't keep track of whether the DD visit was actually
                # executed, so just indicate that we scheduled it
                self.DDVisit = None


    def _schedulePath(self, path, nightNum, filt):
        for visitPair in path:

            """
            # check time to see if we're too ahead or behind schedule
            curTime = self.context.time()
            if abs(curTime - expectedNextVisitTime) > #TODO:

            # check clouds to see if we need to adjust
            areCloudsAround = self.context.areCloudsAround()
            cloudMap = self.context.getCloudMap()
            """

            if visitPair.visit1.isComplete and visitPair.visit2.isComplete:
                # we need this check since a visitPair can be pending
                # even if its first visit has been completed
                # (this happens if weather prevents us from doing
                #  the revisit)
                continue

            # yield an actual visit, not a visitPair, to the telescope
            (ra, dec) = (visitPair.ra, visitPair.dec)
            expTime = Config.WFDExpTime

            # make sure the night won't end before this exposure completes
            twilEnd = sky.twilEnd(Config.surveyStartTime, nightNum)
            if self.context.time() + expTime > twilEnd:
                return

            if not visitPair.visit1.isComplete:
                visitPair.visit1.expTime = expTime
                visitPair.visit1.filter = filt
                yield visitPair.visit1
            else:
                visitPair.visit2.expTime = expTime
                visitPair.visit2.filter = filt
                yield visitPair.visit2

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
            self.tonightsMini.discard(visitPair)

            # keep track of our coverage in each direction
            if self._directionOfDec(visitPair.dec) == NORTH:
                self.NVisitsComplete += 1
            elif self._directionOfDec(visitPair.dec) == SOUTH:
                self.SVisitsComplete += 1
            elif self._directionOfDec(visitPair.dec) == EAST:
                self.EVisitsComplete += 1
            else:
                raise RuntimeError("Completed visit " + str(visit) + \
                                   " is in unknown direction")

    def notifyNightEnd(self):
        if len(self.curNightSlewTimes) == 0:
            return
        if self.nightDirection == NORTH:
            self.NEstAvgSlewTime = np.mean(self.curNightSlewTimes)
        elif self.nightDirection == SOUTHEAST:
            self.SEEstAvgSlewTime = np.mean(self.curNightSlewTimes)
        self.makeupVPs.update(self.tonightsMini)
        self.tonightsMini = None
        filtersequence.maxFChangesNightOver(self.nightDirection)


    def _scheduleRestOfNight(self):
        # TODO questionable function but might need in case of clouds?
        curTime = self.context.time()

    def _getNightRaRange(self, nightNum):
        twilStart = sky.twilStart(Config.surveyStartTime, nightNum)
        twilEnd = sky.twilEnd(Config.surveyStartTime, nightNum)

        raStart = sky.raOfMeridian(twilStart)
        raEnd = sky.raOfMeridian(twilEnd)

        return (raStart, raEnd)

    def _directionOfDec(self, dec):
        if dec > Config.maxDec or dec < Config.minDec:
            raise ValueError("Provided dec of " + str(dec) + " is outside " + \
                             "of the survey area.")
        if dec > self.telescope.latitude + Config.zenithBuffer:
            return NORTH
        elif dec > self.telescope.latitude - Config.zenithBuffer:
            return EAST
        else:
            return SOUTH

    def _getMakeupVPs(self, nightNum, direction):
        # return all the visitPairs in makeupVPs that have
        # yet to be completed 
        raRange = self._getNightRaRange(nightNum)
        if direction == EAST:
            raRange = (raRange[0] + Config.zenithBufferOffset,
                       raRange[1] + Config.zenithBufferOffset)

        visitPairs = set()
        for visitPair in self.makeupVPs:
            if (not visitPair.visit1.isComplete or 
                not visitPair.visit2.isComplete):
                # note that this means a visit returned by this function
                # as pending might already be half-completed.
                # schedule() must handle this case
                if Utils.isRaInRange(visitPair.ra, raRange) and \
                        self._directionOfDec(visitPair.dec) == direction:
                    visitPairs.add(visitPair)
        return visitPairs

    def _getNewMiniSurvey(self, nightNum, direction):
        # called at the beginning of a night when there aren't 
        # enough pending visits left

        (raStart, raEnd) = self._getNightRaRange(nightNum)
        #print "new mini on night", nightNum, "with ra: (", raStart, ",", raEnd, ")"

        newVisitPairs = MiniSurvey.newMiniSurvey(self.telescope, raStart, raEnd, direction)
        return newVisitPairs

    def _generateScans(self, nightNum, visitPairs, direction):
        # partitions the visits into N/S or E/W scans
        # which should each take approximately the same amount
        # of time to execute

        # the scans should be about as wide as the typical 
        # vertical distance between pointings: if they are a
        # different width, the simple algorithm of visiting 
        # pointings in order of increasing/decreasing dec will 
        # lead to longer slew times

        """ here's one way you could go about optimizing:
        this may be problematic though for scans which
        are highly obscured by the moon, since those would
        get super wide

        # the algorithm proceeds iteratively: start with 
        # uniform-width scans and change the widths until the
        # total time (exposure plus slew) is about the same
        # for each scan
        """

        # TODO for now, just partition the visits into 
        # scans of constant width
        # this isn't completely stupid because straying
        # from the meridian by a small amount isn't the 
        # end of the world, and scans of constant width
        # should minimize slew times within the scans 

        raMin, raMax = self._getNightRaRange(nightNum)
        raRange = (raMax - raMin) % (2*np.pi)
        if direction == EAST:
            # easterly scans are offset past the zenith buffer zone
            # in RA by zenithBufferOffset
            ERaMin = raMin + Config.zenithBuffer + Config.zenithBufferOffset
            ERaMax = raMax + Config.zenithBuffer + Config.zenithBufferOffset
            ERaMin %= 2*np.pi
            ERaMax %= 2*np.pi

            validVisitPairs = [v for v in visitPairs
                               if Utils.isRaInRange(v.ra, (ERaMin, ERaMax)) and
                               self._directionOfDec(v.dec) == direction]

            numCols = int(raRange / Config.EScanWidth)
            # there are 2*numCols total scans in a 2xnumCols grid
            """
                       ______________________________________
              __      |  Scan 1  |  Scan 3  |  Scan 5  | ...
             /  \     |          |          |          |
             |  |     --------------------------------------
             \__/     |  Scan 2  |  Scan 4  |  Scan 6  | ...
                      |__________|__________|__________|_____

            """
            dRa = raRange / numCols
            scans = [set() for i in range(numCols*2)]
            for v in validVisitPairs:
                col = int(((v.ra - ERaMin) % (2*np.pi)) / dRa)
                if col == numCols:
                    col -= 1
                scanId = col * 2
                if v.dec < self.telescope.latitude:
                    scanId += 1
                if scanId > len(scans):
                    print scanId, ERaMin, ERaMax, v.ra, v.dec
                scans[scanId].add(v)
            return scans


        if direction == NORTH or direction == SOUTH:
            raMid = raMin + raRange / 2

            validVisitPairs = [v for v in visitPairs 
                               if Utils.isRaInRange(v.ra, (raMin, raMax)) and
                                  self._directionOfDec(v.dec) == direction]

            # the field of view width should be approximately the average 
            # vertical distance between pointings in a scan unless we're
            # very behind in one part of the sky
            # TODO decide how wide to make scans
            scanWidth = 2 * self.telescope.fovWidth
            # mod by 2pi in case raMax < raMin (if the range crosses ra=2pi)
            numScans = raRange / scanWidth
            numScans = int(numScans)
            # make sure we have a multiple of N so we can do N visit scans
            # followed by N revisit scans
            # TODO parametrize N
            N = 2
            numScans = int(N * round(numScans / N))

            # how to divide the sky into scans: cartesian or spherical
            divisionType = "spherical"
            if divisionType == "cartesian":
                sphericalCoords = [[visitPair.ra, visitPair.dec]
                                   for visitPair in validVisitPairs]
                sphericalCoords = np.array(sphericalCoords)
                cartesianCoords = Utils.spherical2Cartesian(sphericalCoords[:,0],
                                                            sphericalCoords[:,1])

                cartesianMid = Utils.spherical2Cartesian(raMid, 0)
                displacementVector = np.cross(cartesianMid, [0,0,-1])
                displacements = np.dot(displacementVector, cartesianCoords.T)

                (minD, maxD) = (displacements.min(), displacements.max())
                normedDisplacements = (displacements - minD) / (maxD - minD)

                # scans are ordered by increasing ra
                scans = [set() for i in range(numScans)]
                properScans = (normedDisplacements * numScans).astype(int)
                # handle the case when normedDisplacements == maxD
                properScans[np.where(properScans == numScans)] = numScans - 1

                for i, visitPair in enumerate(validVisitPairs):
                    scans[properScans[i]].add(visitPair)

            elif divisionType == "spherical":
                scans = [set() for i in range(numScans)]
                for visitPair in validVisitPairs:
                    dRa = (visitPair.ra - raMin) % (2*np.pi)
                    scanId = int(dRa / raRange * numScans)
                    if scanId == numScans:
                        scanId = numScans - 1
                    scans[scanId].add(visitPair)

            return scans
