from __future__ import division
import numpy as np
import Config
from OngoingSurvey import OngoingSurvey
import AstronomicalSky
from Telescope import Telescope
from MiniSurvey import MiniSurvey
from Visit import Visit
from Visit import VisitPair
import Utils

from matplotlib import pyplot as plt
import copy


# add a new mini-survey whenever we dip below this much times
# the expected number of visits needed in a night
NUM_VISITS_PADDING = 1

NORTH = 0
SOUTH = 1
EAST = 2
WEST = 3
SOUTHEAST = 4

class Scheduler:
    def __init__(self, telescope, context):
        self.telescope = telescope
        self.context = context
        # or maybe load from .pkl
        self.ongoingSurvey = set()

        # these estimates get updated at the end of each night
        # (currently estAvgExpTime is not updated)
        self.NEstAvgSlewTime = 7 # seconds
        self.SEEstAvgSlewTime = 7 # seconds
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

        # keep track of which filter we should start out in when looking
        # in the N, S, and E
        self.startFilterIds = {NORTH: 0, SOUTH: 0, EAST: 0}
        self.leftOutFilterIds = {NORTH: 0, SOUTH: 0, EAST: 0}

    def _getFilterIds(self, direction, nScans):
        # note: this method modifies instance state
        # (start filters and left out filters)

        # helper to increment mod len(Telescope.filters)
        def inc(filterId):
            return (filterId + 1) % len(Telescope.filters)

        # decide which filters to use for the night
        leftOut = self.leftOutFilterIds[direction]
        nightsFilterIds = range(0, leftOut) + \
                          range(leftOut + 1, len(Telescope.filters))
        self.leftOutFilterIds[direction] = inc(self.leftOutFilterIds[direction])

        filterId = self.startFilterIds[direction]
        filterIds = []
        for i in range(nScans):
            if filterId not in nightsFilterIds:
                filterId = inc(filterId)
            filterIds.append(filterId)
            if i % 2 == 0:
                filterId = inc(filterId)

        # set up startFilterId for next night
        for i in range(2):
            if self.startFilterIds[direction] not in nightsFilterIds:
                self.startFilterIds[direction] = inc(self.startFilterIds[direction])
            self.startFilterIds[direction] = inc(self.startFilterIds[direction])

        return filterIds

    def schedule(self):
        for nightNum in range(Config.surveyNumNights):
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
                nightDirection = NORTH
            else:
                nightDirection = SOUTHEAST


            # reset the slew times array
            self.curNightSlewTimes = []
            prevAltaz = None
            prevFilter = self.telescope.filters[0]

            # return each visit prescribed by scheduleNight()
            for visit in self._scheduleNight(nightNum, nightDirection):
                radec = np.array([[visit.ra, visit.dec]])
                altaz = AstronomicalSky.radec2altaz(radec, self.context.time())[0]
                if prevAltaz is not None:
                    slewTime = self.telescope.calcSlewTime(prevAltaz, prevFilter,
                                                           altaz, visit.filter)
                    self.curNightSlewTimes.append(slewTime)
                prevAltaz = altaz
                prevFilter = visit.filter
                yield visit
            # if the night isn't over and we've exhausted the visits in 
            # self._scheduleNight, return None until the next night starts
            while AstronomicalSky.nightNum(self.context.time()) == nightNum:
                yield None
            if len(self.curNightSlewTimes) > 0:
                if nightDirection == NORTH:
                    self.NEstAvgSlewTime = np.mean(self.curNightSlewTimes)
                elif nightDirection == SOUTHEAST:
                    self.SEEstAvgSlewTime = np.mean(self.curNightSlewTimes)

    def _scheduleNight(self, nightNum, nightDirection):
        # figure out how many visits we'll probably do tonight
        nightLength = AstronomicalSky.nightLength(nightNum)
        if nightDirection == NORTH:
            t = self.NEstAvgSlewTime + self.estAvgExpTime
        elif nightDirection == SOUTHEAST:
            t = self.SEEstAvgSlewTime + self.estAvgExpTime
        expectedNumVisits = int(nightLength / t)
        # TODO it might be greater than numVisits / 2 if we expect to
        # lose revisits due to weather. May model this at some point
        expectedNumVPs = int(expectedNumVisits / 2)

        # for each direction, figure out if we need to add a new mini
        def topUpVPs(direction):
            # this helper method tops up VPs in direction, then returns
            # the VPs in that direction
            if direction == NORTH:
                expectedVPsInDir = expectedNumVPs * 1
            elif direction == SOUTH:
                expectedVPsInDir = expectedNumVPs * (self.areaInSouth /
                                       (self.areaInSouth + self.areaInEast))
            elif direction == EAST:
                expectedVPsInDir = expectedNumVPs * (self.areaInEast /
                                       (self.areaInSouth + self.areaInEast))

            pendingVPs = self._calculatePendingVisitPairs(nightNum, direction)
            while len(pendingVPs) < expectedVPsInDir * NUM_VISITS_PADDING:
                self._addNewMiniSurvey(nightNum, direction)
                pendingVPs = self._calculatePendingVisitPairs(nightNum, direction)
            return pendingVPs

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
   
        if nightDirection == NORTH:
            pendingVPs = topUpVPs(NORTH)
            # choose which visitPairs we're going to execute tonight
            # TODO could prioritize doing stale visitPairs
            tonightsVPs = np.random.choice(list(pendingVPs),
                                           size=int(expectedNumVPs), 
                                           replace=False)
            scans = self._generateScans(nightNum, tonightsVPs, NORTH)
            # scanningOrder looks liks [0,1,0,1,2,3,2,3,4,5,4,5,...]
            scanningOrder = (np.arange(len(scans))
                               .reshape(int(len(scans)/2), 2)
                               .repeat(2, axis=0)
                               .flatten())
            scans = [scans[i] for i in scanningOrder]

            # scanDirs are the directions to execute each scan in
            # in this case we start out going north, then south, etc
            scanDirs = np.zeros(len(scans))
            scanDirs[0::2] = NORTH
            scanDirs[1::2] = SOUTH

            filterIds = self._getFilterIds(NORTH, len(scans))

        elif nightDirection == SOUTHEAST:
            SPendingVPs = topUpVPs(SOUTH)
            EPendingVPs = topUpVPs(EAST)
            areaInSouthEast = self.areaInSouth + self.areaInEast
            SExpectedNumVPs = expectedNumVPs * self.areaInSouth / areaInSouthEast
            EExpectedNumVPs = expectedNumVPs * self.areaInEast / areaInSouthEast
            STonightsVPs = np.random.choice(list(SPendingVPs),
                                            size=int(SExpectedNumVPs), 
                                            replace=False)
            ETonightsVPs = np.random.choice(list(EPendingVPs),
                                            size=int(EExpectedNumVPs), 
                                            replace=False)

            #print "S tonight", len(STonightsVPs) / self.areaInSouth
            #print "E tonight", len(ETonightsVPs) / self.areaInEast

            SScans = self._generateScans(nightNum, STonightsVPs, SOUTH)
            EScans = self._generateScans(nightNum, ETonightsVPs, EAST)

            SFilterIds = self._getFilterIds(SOUTH, len(SScans)*2)
            EFilterIds = self._getFilterIds(EAST,  len(EScans)*2)

            #print "S scans", len([v for s in SScans for v in s]) / self.areaInSouth
            #print "E scans", len([v for s in EScans for v in s]) / self.areaInEast

            EScansRas = [[v.ra for v in scan] for scan in EScans]
            EMinRas = np.array([min(ras) for ras in EScansRas])
            # gives three numbers: the min RA from scans 1/2, 
            # from scans 3/4, and 5/6
            EMinGroupRa = EMinRas.reshape(3,2).T.min(axis=0)

            # do as many southern pairs as we can before we have to
            # catch the E fields before they hit the zenith avoidance zone
            zenith = np.array([[np.pi/2, 0]])
            raOfZenith = AstronomicalSky.altaz2radec(zenith, self.context.time())[0,0]
            cutoffRa = raOfZenith + Config.zenithBuffer
            timesLeft = (EMinGroupRa - cutoffRa) * (3600*24) / (2*np.pi)

            avgVisitTime = self.estAvgExpTime + self.SEEstAvgSlewTime
            EScanTimes = np.array(map(len, EScans)) * avgVisitTime

            execTimes = EScanTimes[::2] * 2 + EScanTimes[1::2]

            scans = []
            filterIds = []
            scanDirs = []
            cumTime = 0
            j = 0
            for i in range(0, len(SScans), 2):
                numNewVisits = len(SScans[i]) + len(SScans[i+1])
                cumTime += 2 * numNewVisits * avgVisitTime
                if j <= 2 and cumTime > timesLeft[j] - execTimes[j]:
                    scans += [EScans[2*j], EScans[2*j+1], EScans[2*j], EScans[2*j+1]]
                    scanDirs += [WEST, EAST, WEST, EAST]
                    for _ in range(4):
                        filterIds.append(EFilterIds.pop(0))
                    ENumNewVisits = len(EScans[2*j]) + len(EScans[2*j+1])
                    cumTime += 2 * ENumNewVisits * avgVisitTime
                    #cumTime += execTimes[j]
                    j += 1
                scans += [SScans[i], SScans[i+1], SScans[i], SScans[i+1]]
                scanDirs += [SOUTH, NORTH, SOUTH, NORTH]
                for _ in range(4):
                    filterIds.append(SFilterIds.pop(0))
            if j != 3:
                raise RuntimeWarning("j=" + str(j) + "!")

        for scan, scanDir, filterId in zip(scans, scanDirs, filterIds):
            filter = Telescope.filters[filterId]
            if scanDir == NORTH:
                sortedScan = sorted(scan, key=lambda v: v.dec)
            elif scanDir == SOUTH:
                sortedScan = sorted(scan, key=lambda v: -1*v.dec)
            elif scanDir == EAST:
                sortedScan = sorted(scan, key=lambda v: v.ra)
            elif scanDir == WEST:
                sortedScan = sorted(scan, key=lambda v: -1*v.ra)
            else:
                raise RuntimeError("invalid direction " + str(scanDir))

    
            for visit in self._schedulePath(sortedScan, nightNum, filter):
                yield visit


    def _schedulePath(self, path, nightNum, filter):
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
            expTime = AstronomicalSky.getExpTime(ra, dec, "otherstuff")

            # make sure the night won't end before this exposure completes
            nightEnd = AstronomicalSky.nightEnd(nightNum)
            if self.context.time() + expTime > nightEnd:
                return

            if not visitPair.visit1.isComplete:
                visitPair.visit1.expTime = expTime
                visitPair.visit1.filter = filter
                yield visitPair.visit1
            else:
                visitPair.visit2.expTime = expTime
                visitPair.visit2.filter = filter
                yield visitPair.visit2

    def notifyVisitComplete(self, visit, time):
        if not isinstance(visit, Visit):
            raise TypeError("must pass in Visit() instance")

        if visit.isComplete:
            raise RuntimeError("visit was completed twice")

        visit.isComplete = True
        visit.timeOfCompletion = time

        # check if we need to delete the visitPair from the ongoingSurvey
        visitPair = visit.visitPair
        if visitPair.visit1.isComplete and visitPair.visit2.isComplete:
            # take the visit out of our ongoing survey
            self.ongoingSurvey.remove(visitPair)

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

    def _scheduleRestOfNight(self):
        # TODO questionable function but might need in case of clouds?
        curTime = self.context.time()

    def _getNightRaRange(self, nightNum):
        nightStart = AstronomicalSky.nightStart(nightNum)
        nightEnd = AstronomicalSky.nightEnd(nightNum)

        raStart = AstronomicalSky.raOfMeridian(nightStart)
        raEnd = AstronomicalSky.raOfMeridian(nightEnd)

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

    def _calculatePendingVisitPairs(self, nightNum, direction):
        # return all the visitPairs in ongoingSurvey which have
        # yet to be completed 
        raRange = self._getNightRaRange(nightNum)
        if direction == EAST:
            raRange = (raRange[0] + Config.zenithBufferOffset,
                       raRange[1] + Config.zenithBufferOffset)

        visitPairs = set()
        for visitPair in self.ongoingSurvey:
            if (not visitPair.visit1.isComplete or 
                not visitPair.visit2.isComplete):
                # note that this means a visit returned by this function
                # as pending might already be half-completed.
                # schedule() must handle this case
                if Utils.isRaInRange(visitPair.ra, raRange) and \
                        self._directionOfDec(visitPair.dec) == direction:
                    visitPairs.add(visitPair)
        return visitPairs

    def _addNewMiniSurvey(self, nightNum, direction):
        # called at the beginning of a night when there aren't 
        # enough pending visits left

        (raStart, raEnd) = self._getNightRaRange(nightNum)
        #print "new mini on night", nightNum, "with ra: (", raStart, ",", raEnd, ")"

        newVisitPairs = MiniSurvey.newMiniSurvey(self.telescope, raStart, raEnd, direction)
        self.ongoingSurvey.update(newVisitPairs)

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
            raMin += Config.zenithBuffer + Config.zenithBufferOffset
            raMax += Config.zenithBuffer + Config.zenithBufferOffset

            validVisitPairs = [v for v in visitPairs
                               if Utils.isRaInRange(v.ra, (raMin, raMax)) and
                               self._directionOfDec(v.dec) == direction]

            # there are six total scans in a 2x3 grid
            """
                       _________________________________
              __      |  Scan 1  |  Scan 3  |  Scan 5  |
             /  \     |          |          |          |
             |  |     ----------------------------------
             \__/     |  Scan 2  |  Scan 4  |  Scan 6  |
                      |__________|__________|__________|

            """
            dRa = raRange / 3
            ra0 = raMin
            ra1 = (raMin + dRa) % (2*np.pi)
            ra2 = (raMin + 2*dRa) % (2*np.pi)
            ra3 = raMax
            def assignScans(minRa, maxRa, isNorth):
                scan = [v for v in visitPairs 
                        if Utils.isRaInRange(v.ra, (minRa, maxRa)) and
                         ( (v.dec > self.telescope.latitude and isNorth) or
                           (v.dec < self.telescope.latitude and not isNorth) )]
                return set(scan)
            scans = [assignScans(ra0, ra1, True),
                     assignScans(ra0, ra1, False),
                     assignScans(ra1, ra2, True),
                     assignScans(ra1, ra2, False),
                     assignScans(ra2, ra3, True),
                     assignScans(ra2, ra3, False)]
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
