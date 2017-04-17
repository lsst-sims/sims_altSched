from __future__ import division
import numpy as np
import Config
from OngoingSurvey import OngoingSurvey
import AstronomicalSky
from MiniSurvey import MiniSurvey
from Visit import Visit
from Visit import VisitPair
import Telescope
import Utils


# add a new mini-survey whenever we dip below this much times
# the expected number of visits needed in a night
NUM_VISITS_PADDING = 1

NORTH = 0
SOUTH = 1

class Scheduler:
    def __init__(self, context):
        self.context = context
        # or maybe load from .pkl
        self.ongoingSurvey = OngoingSurvey()

        # these estimates get updated at the end of each night
        # (currently estAvgExpTime is not updated)
        self.estAvgSlewTime = 5 # seconds
        self.estAvgExpTime = 30 #seconds

        # keep track of how long we spend slewing each night
        # so we can update our estimate of the average slew time each night
        self.curNightSlewTimes = []

        # int_{telLat+5 deg}^{maxDec} 2\pi r dr, where r=\cos\theta
        # = 2\pi (\sin(maxDec) - \sin(telLat+5 deg)) 
        self.areaInNorth = 2*np.pi*(np.sin(Config.maxDec) - \
                                    np.sin(Telescope.latitude + np.radians(5)))

        # similar for south
        self.areaInSouth = 2*np.pi*(np.sin(Telescope.latitude - np.radians(5)) - \
                                    np.sin(Config.minDec))

        self.numVisitsScheduledInNorth = 0
        self.numVisitsScheduledInSouth = 0

        self.numVisitsPendingInNorth = 0
        self.numVisitsPendingInSouth = 0

    def schedule(self):
        for nightNum in range(Config.surveyNumNights):
            # reset the slew times array
            self.curNightSlewTimes = []
            prevAltaz = None

            # return each visit prescribed by scheduleNight()
            for visit in self._scheduleNight(nightNum):
                radec = np.array([[visit.ra, visit.dec]])
                altaz = AstronomicalSky.radec2altaz(radec, self.context.time())[0]
                if prevAltaz is not None:
                    slewTime = Telescope.calcSlewTime(prevAltaz, altaz)
                    self.curNightSlewTimes.append(slewTime)
                prevAltaz = altaz
                yield visit
            # if the night isn't over and we've exhausted the visits in 
            # self._scheduleNight, return None until the next night starts
            while AstronomicalSky.nightNum(self.context.time()) == nightNum:
                yield None
            self.estAvgSlewTime = np.mean(self.curNightSlewTimes)

    def _scheduleNight(self, nightNum):
        # decide which way to point for the night
        if (self.numVisitsScheduledInNorth / self.areaInNorth < 
                self.numVisitsScheduledInSouth / self.areaInSouth):
            direction = NORTH
        else:
            direction = SOUTH

        # figure out how many visits we'll probably do tonight
        nightLength = AstronomicalSky.nightLength(nightNum)
        t = self.estAvgSlewTime + self.estAvgExpTime
        expectedNumVisits = int(nightLength / t)
        # TODO what if it's not just / 2? <-- why wouldn't it be?
        expectedNumVisitPairs = int(expectedNumVisits / 2)

        # figure out whether we need to add a new mini
        pendingVisitPairs = self._calculatePendingVisitPairs(nightNum, direction)

        while len(pendingVisitPairs) < expectedNumVisitPairs * NUM_VISITS_PADDING:
            self._addNewMiniSurvey(nightNum, direction)
            pendingVisitPairs = self._calculatePendingVisitPairs(nightNum, direction)

        # choose which visitPairs we're going to execute tonight
        # TODO could prioritize doing stale visitPairs
        pendingVisitPairs = list(pendingVisitPairs)
        tonightsVisitPairs = np.random.choice(pendingVisitPairs, 
                                              size=expectedNumVisitPairs, 
                                              replace=False)

        scans = self._generateScans(nightNum, tonightsVisitPairs)
        numScans = len(scans)

        for i in range(0, numScans, 2):
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

            upScan   = sorted(scans[i],   key=lambda visit: -1*visit.dec)
            downScan = sorted(scans[i+1], key=lambda visit:    visit.dec)

            path = upScan + downScan + upScan + downScan

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
                    yield visitPair.visit1
                else:
                    visitPair.visit2.expTime = expTime
                    yield visitPair.visit2

    def notifyVisitComplete(self, visit, time):
        if not isinstance(visit, Visit):
            raise TypeError("must pass in Visit() instance")

        if visit.isComplete:
            raise RuntimeError("visit was completed twice")

        visit.isComplete = True
        visit.timeOfCompletion = time

        # check if we need to notify the minisurvey of a completed visitPair
        visitPair = visit.visitPair
        if visitPair.visit1.isComplete and visitPair.visit2.isComplete:
            # the VisitPair is complete so decrement numIncompleteVisitPairs
            if visitPair.miniSurvey.numIncompleteVisitPairs <= 0:
                raise RuntimeError("negative numIncompleteVisitPairs")
            visitPair.miniSurvey.numIncompleteVisitPairs -= 1
            northDecRange = self._getDecRange(NORTH)
            if northDecRange[0] <= visitPair.dec <= northDecRange[1]:
                self.numVisitsPendingInNorth -= 1
            else:
                self.numVisitsPendingInSouth -= 1

    def _scheduleRestOfNight(self):
        # TODO questionable function but might need in case of clouds?
        curTime = self.context.time()

    def _getNightRaRange(self, nightNum):
        nightStart = AstronomicalSky.nightStart(nightNum)
        nightEnd = AstronomicalSky.nightEnd(nightNum)

        # TODO explain why backwards
        raStart = AstronomicalSky.raOfMeridian(nightStart)
        raEnd = AstronomicalSky.raOfMeridian(nightEnd)

        return (raStart, raEnd)

    def _getDecRange(self, direction):
        if direction == NORTH:
            decRange = (Telescope.latitude + np.radians(5), Config.maxDec)
        else:
            decRange = (Config.minDec, Telescope.latitude - np.radians(5))
        return decRange

    def _calculatePendingVisitPairs(self, nightNum, direction):
        # TODO this method is very slow! Need to run line_profiler 
        # problem is that miniSurveys are never deleted from ongoingSurvey
        # and visitPairs are never removed from miniSurveys

        # return all the visitPairs in ongoingSurvey which have
        # yet to be completed 
        raRange = self._getNightRaRange(nightNum)
        decRange = self._getDecRange(direction)

        visitPairs = set()
        #print "self.ongoingSurvey.miniSurveys length", len(self.ongoingSurvey.miniSurveys)
        for miniSurvey in self.ongoingSurvey.miniSurveys:
            if miniSurvey.numIncompleteVisitPairs > 0:
                # print "miniSurvey.visitPairs length", len(miniSurvey.visitPairs), miniSurvey.numIncompleteVisitPairs
                for visitPair in miniSurvey.visitPairs:
                    if (not visitPair.visit1.isComplete or 
                        not visitPair.visit2.isComplete):
                        # note that this means a visit returned by this function
                        # as pending might already be half-completed.
                        # schedule() must handle this case
                        if Utils.isRaInRange(visitPair.ra, raRange) and \
                           decRange[0] < visitPair.dec < decRange[1]:
                            visitPairs.add(visitPair)
        return visitPairs

    def _addNewMiniSurvey(self, nightNum, direction):
        # called at the beginning of a night when there aren't 
        # enough pending visits left

        (raStart, raEnd) = self._getNightRaRange(nightNum)
        #print "new mini on night", nightNum, "with ra: (", raStart, ",", raEnd, ")"

        minDec, maxDec = self._getDecRange(direction)

        newMini = MiniSurvey.newMiniSurvey(minDec, maxDec, raStart, raEnd)
        numVisitPairs = len(newMini.visitPairs)
        if direction == NORTH:
            self.numVisitsScheduledInNorth += numVisitPairs
            self.numVisitsPendingInNorth += numVisitPairs
        else:
            self.numVisitsScheduledInSouth += numVisitPairs
            self.numVisitsPendingInSouth += numVisitPairs

        self.ongoingSurvey.miniSurveys.append(newMini)

    def _generateScans(self, nightNum, visitPairs):
        # partitions the visits into vertical scans
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

        (raMin, raMax) = self._getNightRaRange(nightNum)
        raMid = raMin + ((raMax - raMin) % (2*np.pi)) / 2

        validVisitPairs = [v for v in visitPairs 
                           if Utils.isRaInRange(v.ra, (raMin, raMax))]

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

        # the field of view width should be approximately the average 
        # vertical distance between pointings in a scan unless we're
        # very behind in one part of the sky
        fovWidth = Telescope.fovWidth
        # TODO decide how wide to make scans
        scanWidth = 2 * fovWidth
        # mod by 2pi in case raMax < raMin (if the range crosses ra=2pi)
        numScans = ((raMax - raMin) % (2*np.pi)) / scanWidth
        numScans = int(numScans)
        # make sure we have a multiple of N so we can do N visit scans followed
        # by N revisit scans
        # TODO parametrize N
        N = 2
        numScans = int(N * round(numScans / N))

        # scans are ordered by increasing ra
        scans = [set() for i in range(numScans)]
        properScans = (normedDisplacements * numScans).astype(int)
        # handle the case when normedDisplacements == maxD
        properScans[np.where(properScans == numScans)] = numScans - 1

        for i, visitPair in enumerate(validVisitPairs):
            scans[properScans[i]].add(visitPair)

        for scan in scans:
            #print len(scan)
            pass
        return scans
