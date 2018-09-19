from __future__ import division
from __future__ import print_function

from config import NORTH, SOUTH, EAST, WEST, SOUTHEAST
import config
import utils
from minis.MiniSurvey import MiniSurvey
from Visit import Visit
from Visit import PROP_WFD, PROP_DD
from minis import filtersequence

from lsst.sims.speedObservatory import sky
from lsst.sims.speedObservatory.utils import unix2mjd

import csv
import numpy as np
import numpy.lib.recfunctions as rf

from DDF import DDScheduler

# initial estimates of slew time (updated after the first night)
_estAvgSlewTimes = {NORTH: 10,
                    SOUTHEAST: 10}

class NightScheduler:
    """ Class that schedules a single night

    The main method is schedule(), though the various notify...() methods
    must be called at the appropriate times as well. See comment in schedule()
    for information on how the NightScheduler schedules.
    """
    def __init__(self, context, telescope, nightNum, direction, makeupVPs):
        self.context = context
        self.telescope = telescope
        self.nightNum = nightNum
        self.direction = direction
        self.makeupVPs = makeupVPs

        # flag that is raised when the dome was just closed
        self.wasDomeJustClosed = False
        self.domeClosedDuration = 0


    def _isDDVisitInRevisitGroup(self, revisitGroup):
        # TODO probably delete this
        # helper to check if a DD visit falls within the min/max
        # RA range of a revisit group
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
        return utils.isRaInRange(self.tonightsDDF['ra'], (minRa, maxRa))


    def schedule(self):
        """ Schedule the night

        This method is a generator that yields Visit instances until it has
        nowhere left to recommend pointing for the night. If self.direction
        is NORTH, it schedules visits in a North/South scanning pattern shown
        in the ascii art comment below. If self.direction is SOUTHEAST, it
        alternates between scheduling visits in the South and visits in the
        zenith dec band. For Southern scans, the pattern is the same as in the
        North. For scans in the zenith dec band (Eastern scans), the scanning
        pattern is similar but goes East/West.

        The NightScheduler chooses which pointings to use for the night by
        generating a new altsched.minis.MiniSurvey and then by drawing from
        self.makeupVPs as necessary if there are not enough visits in the
        new MiniSurvey. It is therefore important to choose a tiling with
        a density that will allow the scheduler to keep up with the meridian
        without drawing too heavily on self.makeupVPs. Note that it is not
        currently supported to observe only a subset of each night's
        MiniSurvey, so if the tiling is too dense, the scheduler will always
        fall behind the meridian.

        Yields
        ------
        The next Visit to be scheduled
        """
        # this will hold the VPs in tonight's mini survey
        self.tonightsMini = None
       
        # figure out how many visits we'll probably do tonight
        nightLength = sky.nightLength(config.surveyStartTime, self.nightNum)

        if config.useDD:
            # figure out how much time we'll spend doing DD tonight
            # TODO bad accessing private thing -- just make a getter...
            leftOutFilter = filtersequence._leftOutFilters[self.direction]
            ddScheduler = DDScheduler(self.nightNum, self.direction, leftOutFilter)
            ddTime = ddScheduler.telTimeTonight()
            ddDesiredStartTime = ddScheduler.desiredStartTime()
        else:
            ddTime = 0


        # the period of exposure taking is the sum of the exposure time,
        # the visit overhead time (shutter/intermediate readout),
        # and the slew time
        expectedWfdVisitTime = self._expectedWfdVisitTime()
        expectedNumVisits = int((nightLength - ddTime) / expectedWfdVisitTime)

        """
        if self.tonightsDDF is None:
            # we won't do a DD visit tonight
            expectedNumVisits = int(nightLength / expectedWfdVisitTime)
        else:
            # we do plan to do a DD visit tonight, so the amount of time
            # available for WFD exposures is nightLength - DDexpTime
            # TODO this doesn't include slew time to/from the DD field
            expTime = self.tonightsDDF.sequences[self.tonightsDDF.currentSequence].totExpTime

            expectedNumVisits = int((nightLength - ddTime) / expectedWfdVisitTime)
        """

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


        # a revisitGroup is a pair of scans with information
        # about how to execute them
        if self.direction == NORTH:
            revisitGroups = self._getNorthRevisitGroups(expectedNumVPs)
        elif self.direction == SOUTHEAST:
            revisitGroups = self._getSoutheastRevisitGroups(expectedNumVPs)

        # now schedule each revisit group
        rgIdx = 0
        # this is used in case the dome closes during a revisitgroup to track
        # how much observing time we would have accomplished during the rest
        # of that revisitgroup
        numVisitsKilled = 0

        isDDComplete = True if ddTime == 0 else False
        while rgIdx < len(revisitGroups) or not isDDComplete:
            assert(rgIdx <= len(revisitGroups))

            # first, check if we should do DD now
            if config.useDD:
                # flag to catch the case that the optimal DD time is at the
                # very end of the night
                onlyDDLeft = rgIdx == len(revisitGroups)

                # flag to check if we should execute the DD before doing the next
                # revisit group
                isTimeForDD = False
                if rgIdx < len(revisitGroups):
                    nextRg = revisitGroups[rgIdx]
                    nVisits = 2 * (len(nextRg["scan1"]) + len(nextRg["scan2"]))
                    timeForNextBlock = nVisits * expectedWfdVisitTime
                    curTime = self.context.time()
                    # we should do the DD now if executing the next revisit group
                    # would make us late for the DD
                    isTimeForDD = curTime + timeForNextBlock / 2 >= ddDesiredStartTime
                if (not isDDComplete) and (onlyDDLeft or isTimeForDD):
                    # figure out which filters we should suggest that the DD
                    # observations start and end in
                    if rgIdx == 0:
                        startFilter = None
                    else:
                        startFilter = revisitGroups[rgIdx - 1]["filter4"]
                    if rgIdx == len(revisitGroups):
                        endFilter = None
                    else:
                        endFilter = revisitGroups[rgIdx]["filter1"]
                    # TODO start/end filter don't make sense as is since
                    # they'll always be the same (for maxfchanges). Should
                    # instead modify WFD so it does e.g. g-r-i-z instead of
                    # g-r-r-i surrounding the DD observation. Or just switch
                    # to 2 filters per night, which might make things easier
                    # (then you'd do, e.g. g-r-g-r instead of g-r-r-g around
                    # a dd observation)
                    for visit in ddScheduler.schedule(startFilter, endFilter):
                        yield visit
                        if self.wasDomeJustClosed:
                            # assumes we got through none of the DD before the dome
                            # closed, which won't usually be true but doesn't
                            # really matter and saves the hassle of keeping track
                            # of how much did complete

                            # also assumes that the dome closed for at least the
                            # length of the DD visit

                            # note that we don't reset self.wasDomeJustClosed
                            # because if the dome was closed for a while, we should
                            # also kill some revisit groups (done below)
                            self.domeClosedDuration -= ddTime
                    isDDComplete = True
                if rgIdx == len(revisitGroups):
                    # we just had the DD left; no more WFD revisit groups
                    break

            rg = revisitGroups[rgIdx]
            rgIdx += 1

            # perform scan1 and scan2 twice in a row but use
            # potentially-different filters the second time through
            for i, scan, scanDir, filt in [
                    (0, rg["scan1"], rg["scanDir1"], rg["filter1"]),
                    (1, rg["scan2"], rg["scanDir2"], rg["filter2"]),
                    (2, rg["scan1"], rg["scanDir1"], rg["filter3"]),
                    (3, rg["scan2"], rg["scanDir2"], rg["filter4"])]:

                if scanDir == NORTH:
                    sortedScan = sorted(scan, key=lambda v: v.dec)
                elif scanDir == SOUTH:
                    sortedScan = sorted(scan, key=lambda v: -1*v.dec)
                elif scanDir == EAST:
                    (raMin, raMax) = self._getTonightsRaRange()
                    # subtract raMin so we never cross 2\pi (since a scan
                    # is never more than 2\pi radians long in RA)
                    sortedScan = sorted(scan, key=lambda v: (v.ra - raMin) % (2*np.pi))
                elif scanDir == WEST:
                    (raMin, raMax) = self._getTonightsRaRange()
                    sortedScan = sorted(scan, key=lambda v: -1 * (v.ra - raMin) % (2*np.pi))
                else:
                    raise RuntimeError("invalid direction " + str(scanDir))

                for visit in self._schedulePath(sortedScan, filt):
                    # _schedulePath returns before yielding a stale visit
                    # due to a closed dome. So if it returns a visit, we
                    # know the dome wasn't just closed
                    yield visit

                if self.wasDomeJustClosed:
                    # kill the rest of this revisit group if the dome was just
                    # closed
                    # TODO doesn't consider time killed in the current scan
                    if i == 0:
                        numVisitsKilled += len(rg["scan1"]) + 2 * len(rg["scan2"])
                    elif i == 1:
                        numVisitsKilled += len(rg["scan1"]) + len(rg["scan2"])
                    elif i == 2:
                        numVisitsKilled += len(rg["scan2"])

                    break

            if self.wasDomeJustClosed:
                avgVisitTime = self._expectedWfdVisitTime()
                timeKilled = 0
                while rgIdx < len(revisitGroups) - 1:
                    rg = revisitGroups[rgIdx]
                    # this assumes we completed no visits from the RG during
                    # which the dome closed (which won't usually be true, but
                    # doesn't matter much and saves us the hassle of calculating
                    # how many did complete successfully)
                    numVisitsKilled += (len(rg["scan1"]) + len(rg["scan2"])) * 2

                    # decide whether or not to kill one more revisitgroup
                    stopNow = np.abs(timeKilled - self.domeClosedDuration)
                    oneMore = np.abs(numVisitsKilled * avgVisitTime -
                                     self.domeClosedDuration)
                    if stopNow < oneMore:
                        break
                    else:
                        timeKilled = numVisitsKilled * avgVisitTime
                        rgIdx += 1

                # we finished handling the dome closed interrupt, so lower the
                # flag
                self.wasDomeJustClosed = False
        # end loop over revisit groups
        if config.useDD:
            assert(isDDComplete)

    def scheduleDD(self):
        """ Schedule DD visits """

        ra = self.tonightsDDF.ra
        dec = self.tonightsDDF.dec
        currentSequence = self.tonightsDDF.currentSequence
        for visit in self.tonightsDDF.sequences[currentSequence].schedule(ra, dec):
            yield visit

        self.tonightsDDF.prevNight = self.nightNum
        self.tonightsDDF.last_sequence = currentSequence
        
    def _schedulePath(self, path, filt):
        """ Schedules a list of pointings

        This method only schedules one visit of each passed-in visitPair.
        """
        for visitPair in path:
            # check time to see if the dome closed since the last
            # visit. If so, the rest of this path is stale
            # TODO "path" is too generic since a path could be arbitrarily long
            if self.wasDomeJustClosed:
                return

            if visitPair.visit1.isComplete and visitPair.visit2.isComplete:
                # we need this check since a visitPair can be pending
                # even if its first visit has been completed
                # (this happens if weather prevents us from doing
                #  the revisit)
                continue

            # yield an actual visit, not a visitPair, to the telescope
            (ra, dec) = (visitPair.ra, visitPair.dec)
            expTime = config.WFDExpTime

            if not visitPair.visit1.isComplete:
                visitPair.visit1.expTime = expTime
                visitPair.visit1.filter = filt
                yield visitPair.visit1
            else:
                visitPair.visit2.expTime = expTime
                visitPair.visit2.filter = filt
                yield visitPair.visit2


    def notifyDomeClosed(self, durationClosed):
        """ Method used to notify us that the dome just closed

        This method should be called every time the dome is closed
        or else the NightScheduler won't know to skip ahead. Note that the
        schedule() method could check the current time after every visit
        to make sure that not too much time has elapsed. But putting the
        semantics of dome closing into an interrupt routine like this
        makes explicit the control system/scheduler communication about
        the dome closing that would be implicit if we just checked the time.

        Parameters
        ----------
        durationClosed : float
            The number of seconds that the dome was closed for
        """
        # wasDomeJustClosed is an interrupt flag for schedule
        self.wasDomeJustClosed = True
        self.domeClosedDuration = durationClosed

    def notifyVisitPairComplete(self, visitPair):
        """ Method used to notify us that a VisitPair was carried out

        This should be called every time a VisitPair is completed (i.e.
        both constituent Visits were carried out).

        Parameters
        ----------
        visitPair : altSched.VisitPair
            The VisitPair that was just completed
        """
        self.tonightsMini.discard(visitPair)

    def notifyNightEnd(self, avgSlewTime):
        """ Method used to notify us that the night is over

        This method updates the average slew time and manages
        the filter bookkeeping

        Parameters
        ----------
        avgSlewTime : float
            The average slew time during the night that just ended (seconds)

        Returns
        -------
        A set of VisitPairs that were part of tonight's tiling but did not
        get fully executed.
        """
        # update the _estAvgSlewTimes for tonight's direction
        _estAvgSlewTimes[self.direction] = avgSlewTime
        filtersequence.maxFChangesNightOver(self.direction)
        return self.tonightsMini

    def _getMakeupVPs(self, direction):
        # return all the visitPairs in makeupVPs that have
        # yet to be completed
        raRange = self._getTonightsRaRange()
        if direction == EAST:
            raRange = (raRange[0] + config.zenithBufferOffset,
                       raRange[1] + config.zenithBufferOffset)

        visitPairs = set()
        for visitPair in self.makeupVPs:
            if (not visitPair.visit1.isComplete or
                not visitPair.visit2.isComplete):
                # note that this means a visit returned by this function
                # as pending might already be half-completed.
                # schedule() must handle this case
                if utils.isRaInRange(visitPair.ra, raRange) and \
                        utils.directionOfDec(visitPair.dec) == direction:
                    visitPairs.add(visitPair)
        return visitPairs


    def _getFilters(self, nScans):
        return filtersequence.maxFChanges(self.nightNum, self.direction, nScans)

    def _expectedWfdVisitTime(self):
        """ Returns the expected average visit time including slew """
        return (_estAvgSlewTimes[self.direction] +
                config.WFDExpTime +
                config.visitOverheadTime)

    def _getTonightsDDF(self):
        # TODO probably no longer needed
        """ Decides which DD field to observe for tonight """

        nightStart = sky.nightStart(config.surveyStartTime, self.nightNum)
        nightEnd = sky.nightEnd(config.surveyStartTime, self.nightNum)
        minRa = sky.raOfMeridian(nightStart)
        maxRa = sky.raOfMeridian(nightEnd)

        # loop through all DD fields to see which one we can observe tonight

        visibleDDFs = []

        r = []
    
        for ddf in self.DDFs:
            ra = ddf.ra
            dec = ddf.dec
            ddf.currentSequence = ddf.last_sequence + 1
            if ddf.currentSequence >= ddf.nseq:
                ddf.currentSequence = 0
            exptime = ddf.sequences[ddf.currentSequence].totExpTime

            raBuffer = exptime / (24*3600) * 2*np.pi
            DDDirection = utils.directionOfDec(dec)

            if DDDirection == NORTH:
                inRaRange = utils.isRaInRange(ra, (minRa, maxRa - raBuffer))
            elif DDDirection == SOUTH:
                inRaRange = utils.isRaInRange(ra, (minRa, maxRa - raBuffer))
                DDDirection = SOUTHEAST
            elif DDDirection == EAST:
                # DD visits in the East need to be in the offset
                # allowable RA range
                raRange = (minRa + config.zenithBufferOffset + raBuffer,
                           maxRa + config.zenithBufferOffset - raBuffer)
                inRaRange = utils.isRaInRange(ra, raRange)
                DDDirection = SOUTHEAST

            if not inRaRange:
                continue

            # do not take into account fields that were observed
            # recently (2 days ago)

            if ddf.prevNight == -1 or self.nightNum - ddf.prevNight > 2:
                th_min, th_max, ha, time_obs = \
                        utils.getTimeObs(ra, dec, exptime, self.telescope,
                                         nightStart, nightEnd)
                ddf.time_obs = time_obs
                ddf.ha = ha
                visibleDDFs.append(ddf)

        if len(visibleDDFs) > 0:
            # among the selected, if there are fields not observed (yet) ie
            # prevNight = -1 : should be priority
            list_priority=[]
            for ddf in visibleDDFs:
                #print('selected',ddf)
                if ddf.prevNight == -1:
                    list_priority.append(ddf)

            if len(list_priority) > 0:
                visibleDDFs = list_priority

            if len(visibleDDFs) > 1 :
                r = []
                for i, ddf in enumerate(visibleDDFs):
                    r.append((i, ddf.prevNight, ddf.ha))
                tab = np.rec.fromrecords(r, names=['index', 'prevNight', 'ha'])
                """
                idv = tab['prevNight'] == np.min(tab['prevNight'])
                sel_visit=tab[idv]
                if len(sel_visit) > 1:
                """
                min_ha = np.argmin(tab['ha'])
                return visibleDDFs[np.asscalar(tab[min_ha]['index'])]
                """
                else:
                    return visibleDDFs[np.asscalar(sel_visit['index'])]
                """
            else:
                return visibleDDFs[0]
        else:
            return None

    def _getTonightsRaRange(self):
        twilStart = sky.twilStart(config.surveyStartTime, self.nightNum)
        twilEnd = sky.twilEnd(config.surveyStartTime, self.nightNum)

        raStart = sky.raOfMeridian(twilStart)
        raEnd = sky.raOfMeridian(twilEnd)

        return (raStart, raEnd)

    def _getNewMiniSurvey(self, direction):
        # called at the beginning of a night when there aren't
        # enough pending visits left

        (raStart, raEnd) = self._getTonightsRaRange()
        #print("new mini on night", nightNum, "with ra: (", raStart, ",", raEnd, ")")

        newVisitPairs = MiniSurvey.newMiniSurvey(self.telescope, raStart,
                                                 raEnd, direction)
        return newVisitPairs

    def _generateScans(self, visitPairs, direction):
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

        raMin, raMax = self._getTonightsRaRange()
        raRange = (raMax - raMin) % (2*np.pi)
        if direction == EAST:
            # easterly scans are offset past the zenith buffer zone
            # in RA by zenithBufferOffset
            ERaMin = raMin + config.zenithBuffer + config.zenithBufferOffset
            ERaMax = raMax + config.zenithBuffer + config.zenithBufferOffset
            ERaMin %= 2*np.pi
            ERaMax %= 2*np.pi

            validVisitPairs = [v for v in visitPairs
                               if utils.isRaInRange(v.ra, (ERaMin, ERaMax)) and
                               utils.directionOfDec(v.dec) == direction]

            numCols = int(raRange / config.EScanWidth)
            # there are 2*numCols total scans in a 2xnumCols grid
            """
            zenith     ______________________________________
              __      |  Scan 1  |  Scan 3  |  Scan 5  | ...
             /  \     |          |          |          |
             |  |     --------------------------------------
             \__/     |  Scan 2  |  Scan 4  |  Scan 6  | ...
            avoidance |__________|__________|__________|_____

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
                    print(scanId, ERaMin, ERaMax, v.ra, v.dec)
                scans[scanId].add(v)
            return scans


        if direction == NORTH or direction == SOUTH:
            raMid = raMin + raRange / 2

            validVisitPairs = [v for v in visitPairs
                               if utils.isRaInRange(v.ra, (raMin, raMax)) and
                                  utils.directionOfDec(v.dec) == direction]

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
                cartesianCoords = utils.spherical2Cartesian(sphericalCoords[:,0],
                                                            sphericalCoords[:,1])

                cartesianMid = utils.spherical2Cartesian(raMid, 0)
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

    def _getNorthRevisitGroups(self, requestedNumVPs):
        makeupVPs = self._getMakeupVPs(NORTH)
        if len(makeupVPs) >= requestedNumVPs:
            # tonight will be a purely makeup night
            self.tonightsMini = set()
            tonightsVPs = np.random.choice(list(makeupVPs),
                                           size=int(requestedNumVPs),
                                           replace=False)
        else:
            # we start with one new minisurvey for tonight and then
            # fill in the remaining needed visits with makeup VPs
            self.tonightsMini = self._getNewMiniSurvey(NORTH)

            # add to makeupVPs if we still don't have enough visit pairs
            while(len(self.tonightsMini) + len(makeupVPs) <
                    requestedNumVPs):
                self.makeupVPs.update(self._getNewMiniSurvey(NORTH))
                makeupVPs = self._getMakeupVPs(NORTH)
            # add some VPs from the makeups to tonightsVPs
            tonightsVPs = list(self.tonightsMini)
            if len(self.tonightsMini) < int(requestedNumVPs):
                tonightsVPs.extend(
                    np.random.choice(list(makeupVPs),
                                     size=int(requestedNumVPs) - len(self.tonightsMini),
                                     replace=False)
                )

        # choose which visitPairs we're going to execute tonight
        # TODO could prioritize doing stale visitPairs
        scans = self._generateScans(tonightsVPs, NORTH)
        filters = self._getFilters(len(scans)*2)

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

        return revisitGroups

    def _getSoutheastRevisitGroups(self, requestedNumVPs):
        SExpectedNumVPs = (requestedNumVPs * utils.areaInDir(SOUTH) /
                           utils.areaInDir(SOUTHEAST))
        EExpectedNumVPs = (requestedNumVPs * utils.areaInDir(EAST) /
                           utils.areaInDir(SOUTHEAST))

        SMakeupVPs = self._getMakeupVPs(SOUTH)
        EMakeupVPs = self._getMakeupVPs(EAST)
        if(len(SMakeupVPs) >= SExpectedNumVPs and
           len(EMakeupVPs) >= EExpectedNumVPs):
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
            SMini = self._getNewMiniSurvey(SOUTH)
            EMini = self._getNewMiniSurvey(EAST)
            self.tonightsMini = SMini | EMini

            # if needed, add minis to S and E makeupVPs
            while(len(SMini) + len(SMakeupVPs) <
                    SExpectedNumVPs):
                self.makeupVPs.update(self._getNewMiniSurvey(SOUTH))
                SMakeupVPs = self._getMakeupVPs(SOUTH)

            while(len(EMini) + len(EMakeupVPs) <
                    EExpectedNumVPs):
                self.makeupVPs.update(self._getNewMiniSurvey(EAST))
                EMakeupVPs = self._getMakeupVPs(EAST)

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

        SScans = self._generateScans(STonightsVPs, SOUTH)
        EScans = self._generateScans(ETonightsVPs, EAST)

        # TODO confusing use of the word scan to mean both a single vertical
        # strip of visits (that gets executed twice) and a single execution
        # of a verticle strip.
        nSScans = len(SScans)
        nScans = (len(SScans) + len(EScans)) * 2
        filters = self._getFilters(nScans)

        # do as many southern pairs as we can before we have to
        # catch the E fields before they hit the zenith avoidance zone

        raOfZenith = self._getTonightsRaRange()[0]
        cutoffRa = raOfZenith + config.zenithBuffer

        EScansRas = [np.array([v.ra for v in scan]) for scan in EScans]
        EMinRas = np.zeros(len(EScansRas))
        for i, ras in enumerate(EScansRas):
            if len(ras) == 0:
                # TODO
                #print("len(ras)=0 :(. nightNum =", nightNum)
                #print("nightLen", sky.nightLength(config.surveyStartTime, nightNum) / 3600)
                #print("EScans len", [len(scan) for scan in EScans])
                #print("SScans len", [len(scan) for scan in SScans])
                #print("avgVisitTime", self.estAvgExpTime + self.SEEstAvgSlewTime)
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

        avgVisitTime = (config.WFDExpTime + _estAvgSlewTimes[SOUTHEAST] +
                        config.visitOverheadTime)
        EScanTimes = np.array(list(map(len, EScans))) * avgVisitTime

        # this is the time available between when we start a pair of East
        # scans and when we have to be done with those scans. It is NOT
        # the amount of time it would take to complete each pair of East
        # scans (that would be 2 * (EScanTimes[::2] + EScanTimes[1::2]))
        # but rather the time to do the South scan once and the North scan
        # twice, since at that point we're headed back East again so are no
        # longer worried about hitting the zenith avoid zone
        execTimes = 2 * EScanTimes[::2] + EScanTimes[1::2]

        def printDebug():
            print("nightNum", self.nightNum)
            print("numECols", numECols)
            print("EMinGroupRa", EMinGroupRa)
            print("raOfZenith", raOfZenith)
            print("cutoffRa", cutoffRa)
            print("timesLeft", timesLeft, timesLeft / 3600)
            print("avgVisitTime", avgVisitTime)
            print("EScanTimes", EScanTimes, EScanTimes / 3600, np.cumsum(EScanTimes/3600))
            print("execTimes", execTimes)
            SScanTimes = np.array(map(len, SScans)) * avgVisitTime
            print("SScanTimes", SScanTimes, SScanTimes / 3600, np.cumsum(SScanTimes/3600)
)
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
            if config.useDD:
                cumTime += config.DDTimePerNight / nSScans * 2
            paddedTime = cumTime + 0.1 * numNewVisits * avgVisitTime

            # add an East scan if we've run out of time
            # before the scan hits the zenith avoid zone
            if j < numECols and paddedTime > timesLeft[j] - execTimes[j]:
                ERevisitGroup = newRevisitGroup(EAST)
                revisitGroups.append(ERevisitGroup)
                ENumNewVisits = (len(ERevisitGroup["scan1"]) +
                                 len(ERevisitGroup["scan2"]))
                cumTime += 2 * ENumNewVisits * avgVisitTime

                j += 1
            # now add the South scan. Note this assumes that we never need
            # to add two East scans in a row in order to keep out of the
            # zenith avoidance region, which seems like a reasonable assumption
            # TODO but might be false in the case of extreme weather causing
            # a buildup in one southern scan
            SRevisitGroup = newRevisitGroup(SOUTH)

            """
            # add to cumTime if we know we'll schedule the DD visit after
            # this revisit group
            if (self.tonightsDDF is not None and
                    self._isDDVisitInRevisitGroup(SRevisitGroup)):
                # TODO this ignores the slew time to/from the DD field
                cumTime += self.tonightsDDF['totExpTime']
            """
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

        return revisitGroups

