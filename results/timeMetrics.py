from __future__ import division

from lsst.sims.maf.metrics import BaseMetric
import numpy as np
import time
from scipy.optimize import minimize_scalar
from collections import Counter

from lsst.sims.speedObservatory import utils
from lsst.sims.speedObservatory import sky
from lsst.sims.speedObservatory import Telescope

class SNLotsMetric(BaseMetric):
    def __init__(self, nightCol="night", filterCol="filter",
                 m5Col="fiveSigmaDepth", metricName="SNLots", **kwargs):
        super(SNLotsMetric, self).__init__(col=[nightCol, filterCol, m5Col],
                                           units="Fraction", metricName=metricName,
                                           **kwargs)
        self.nightCol = nightCol
        self.filterCol = filterCol
        self.m5Col = m5Col

        self.durationDays = 45
        self.peakDays = 15
        self.riseSlope = -2 / self.peakDays
        self.declineSlope = 1.4 / 30
        self.surveyNYears = 10
        self.peaks = {"g": 23.6, "r": 22.6, "i": 22.7, "z": 22.7}
        self.nFiltersPrePeak = 3
        self.nObsPerFilterPrePeak = 1
        self.nFiltersPostPeak = 4
        self.nObsPerFilterPostPeak = 2


    def run(self, dataSlice, slicePoint=None):
        # Fraction of z=0.5 type Ia SN that are detected in at least 3 filters
        # before peak and twice in all filters post-peak (considering only
        # the griz filters)

        # assume SN go off back to back from night 0 to night
        # self.surveyNYears * 365

        nEvents = int(np.ceil(self.surveyNYears * 365 / self.durationDays))

        # loop over all events, checking if each one is detected
        # enough times to count
        nQualifyingEvents = 0
        for i in range(nEvents):
            eventStart = i * self.durationDays
            eventEnd = eventStart + self.durationDays
            thisEvent = dataSlice[np.where(
                            (eventStart <= dataSlice[self.nightCol]) &
                            (dataSlice[self.nightCol] < eventEnd)
                        )]
            thisEvent[self.nightCol] -= eventStart

            # get detection limits in each band
            detectionLimit = np.zeros(thisEvent[self.nightCol].size)
            for f in self.peaks.keys():
                curF = np.where(thisEvent[self.filterCol] == f)
                detectionLimit[curF] = self.lightCurve(thisEvent[self.nightCol][curF], f)
            assert((detectionLimit > 0).all())

            # figure out which observations detected the thing
            detected = np.where(thisEvent[self.m5Col] >= detectionLimit)

            # count how many were pre/post-peak
            prePeak = np.where(thisEvent[self.nightCol] <= self.peakDays)
            postPeak = np.where(thisEvent[self.nightCol] > self.peakDays)

            prePeakNDetections = Counter(thisEvent[self.filterCol][prePeak])
            postPeakNDetections = Counter(thisEvent[self.filterCol][postPeak])
            assert("u" not in prePeakNDetections and "y" not in prePeakNDetections)
            
            nQualifyingPrePeakFilters = 0
            for f in prePeakNDetections:
                if prePeakNDetections[f] >= self.nObsPerFilterPrePeak:
                    nQualifyingPrePeakFilters += 1

            nQualifyingPostPeakFilters = 0
            for f in postPeakNDetections:
                if postPeakNDetections[f] >= self.nObsPerFilterPostPeak:
                    nQualifyingPostPeakFilters += 1
            if nQualifyingPrePeakFilters >= self.nFiltersPrePeak and \
               nQualifyingPostPeakFilters >= self.nFiltersPostPeak:
                   nQualifyingEvents += 1
        

        return nQualifyingEvents / nEvents

    def lightCurve(self, time, f):
        """
        Calculate the magnitude of the object at each time, in each filter.

        Parameters
        ----------
        time : numpy.ndarray
            The times of the observations (days since the start of the event).
        f: str
            The filter of the observations.

        Returns
        -------
        numpy.ndarray
            The magnitudes of the object at each time, in each filter.
        """

        lcMags = np.ones(time.size, dtype=float) * self.peaks[f]

        rise = np.where(time <= self.peakDays)
        lcMags[rise] += self.riseSlope * (time[rise] - self.peakDays)

        decline = np.where(time > self.peakDays)
        lcMags[decline] += self.declineSlope * (time[decline] - self.peakDays)

        return lcMags


class LonelinessMetric(BaseMetric):
    def __init__(self, dTCutoff, mjdCol="observationStartMJD", **kwargs):
        super(LonelinessMetric, self).__init__(col=[mjdCol], 
                                               units="Fraction", **kwargs)
        self.mjdCol=mjdCol
        self.dTCutoff = dTCutoff

    def run(self, dataSlice, slicePoint=None):
        if dataSlice.size < 2:
            return self.badval

        times = np.sort(dataSlice[self.mjdCol])
        dTs = np.diff(times)

        # dTCuttof needs to be in units of days
        daysCutoff = self.dTCutoff / (3600*24)
        
        # pad dTs with a value over the cutoff (since the first and last
        # visits have no visits before and after, respectively)
        dTs = np.pad(dTs, 1, 'constant', constant_values=daysCutoff+1)
        nLonelies = ((dTs[:-1] > daysCutoff) & (dTs[1:] > daysCutoff)).sum()
        return nLonelies / len(times)

class IntranightColorMetric(BaseMetric):
    def __init__(self, mjdCol="observationStartMJD", nightCol="night",
                 filterCol="filter", **kwargs):
        super(IntranightColorMetric, self).\
              __init__(col=[nightCol, mjdCol, filterCol],
                       units="Fraction", **kwargs)
        self.mjdCol = mjdCol
        self.nightCol = nightCol
        self.filterCol = filterCol

    def run(self, dataSlice, slicePoint=None):
        """ Considering only nights when this point was visited at least twice,
        this metric returns the fraction of nights when a color was obtained
        at this point (i.e. visits in 2+ filters)"""

        if dataSlice.size < 2:
            return self.badval

        order = np.argsort(dataSlice[self.mjdCol])
        nights = dataSlice[self.nightCol][order]
        mjds = dataSlice[self.mjdCol][order]
        filters = dataSlice[self.filterCol][order]

        # number of nights that have at least 2 visits
        numCountedNights = 0

        # number of nights where a color was obtained
        numNightsWithColor = 0

        curNight = 0
        nVisitsTonight = 0
        filtersTonight = set()
        for night, mjd, filt in zip(nights, mjds, filters):
            if night != curNight:
                if nVisitsTonight >= 2:
                    numCountedNights += 1
                    if len(filtersTonight) >= 2:
                        numNightsWithColor += 1
                curNight = night
                nVisitsTonight = 0
                filtersTonight = set()
            nVisitsTonight += 1
            filtersTonight.add(filt)

        if numCountedNights == 0:
            return self.badval

        return numNightsWithColor / numCountedNights

class UnobservedDurationMetric(BaseMetric):
    def __init__(self, cutoffTime, surveyStartTime, surveyEndTime,
                 mjdCol="observationStartMJD", nightCol="night", **kwargs):
        super(UnobservedDurationMetric, self).\
                __init__(col=[nightCol, mjdCol],
                         units="Fraction of survey duration", **kwargs)
        # all are units of seconds
        self.cutoffTime = cutoffTime
        self.surveyStartTime = surveyStartTime
        self.surveyEndTime = surveyEndTime

        self.mjdCol = mjdCol
        self.nightCol = nightCol
    
    def run(self, dataSlice, slicePoint=None):
        """ Calculates the amount of time that this point spends
        unobserved for more than self.cutoffTime as a fraction of the
        complete survey duration. Note this includes time during the day."""

        # convert to units of days
        startTime = utils.unix2mjd(self.surveyStartTime)
        endTime = utils.unix2mjd(self.surveyEndTime)
        cutoffTime = self.cutoffTime / (3600*24)

        runLengthYears = round((endTime - startTime) / 365)

        order = np.argsort(dataSlice[self.mjdCol])
        nights = dataSlice[self.nightCol][order]
        times = dataSlice[self.mjdCol][order]

        # find the night for which it is true that at midnight,
        # this field is closest to opposite the meridian

        # calculating the alt at midnight for every night is slow
        # calculating only once a fortnight is fine since fields
        # always go unobserved for at least two weeks after they
        # set for the year
        oneYear = np.arange(0, 365, 14)
        nightStarts = np.array([sky.nightStart(self.surveyStartTime, night)
                                for night in oneYear])
        nightEnds = np.array([sky.nightEnd(self.surveyEndTime, night)
                              for night in oneYear])
        midnights = (nightEnds + nightStarts) / 2
        # setting dec to Telescope.latitude makes the field pass through nadir

        alts = [sky.radec2altaz(slicePoint["ra"], Telescope.latitude, midnight)[0]
                for midnight in midnights]

        alts = np.array(alts)
        coverageGapCenterNight = oneYear[np.argmin(alts)]

        # calculate how much time the field spends unobserved
        firstObsTime = times[0]
        lastObsTime = times[-1]

        assert(firstObsTime >= startTime)
        assert(lastObsTime  <= endTime)

        dTs = np.diff(times)

        # take out dT values that correspond to times when the field
        # had set for the year
        # loop through every year, finding the first night after the
        # night at the center of the yearly gap. The dT to kill is at
        # that index minus 1

        # keep track of how much time the field spends invisible
        # so we don't include that time in the denominator
        # of the returned fraction

        # note that a field becomes "invisible" precicely after the last
        # observation of that field takes place for the year, and similarly
        # for becoming visible
        invisibleTime = 0
        includeEnds = True
        maxNight = np.ceil(max(nights) /  365) * 365
        while coverageGapCenterNight < maxNight:
            yearlyGapIdx = np.searchsorted(nights, coverageGapCenterNight) - 1
            if yearlyGapIdx >= 0 and yearlyGapIdx < len(nights) - 1:
                # if yearlyGapIdx is -1, that means no visits were carried
                # out this year, so there's no dT to kill

                # if yearlyGapIdx = len(nights) - 1, then searchsorted
                # returned the last night before the gap instead of the
                # first night after the gap
                invisibleTime += dTs[yearlyGapIdx]
                dTs[yearlyGapIdx] = 0
            else:
                # in this case, we should assume the field was
                # invisible until the first observation and* invisible
                # between the last observation and the end time
                # (this will only be wrong for pixels that aren't
                # observed in the first year ever)

                # *Note: this assumes the survey lasts an integer
                # number of years so that whatever is invisible at the
                # beginning is also invisible at the end
                includeEnds = False
                invisibleTime += firstObsTime - startTime
                invisibleTime += endTime - lastObsTime

            coverageGapCenterNight += 365

        if endTime - startTime - invisibleTime == 0:
            # this could happen if there are very few observations
            return self.badval

        # at this point, includeEnds will only be True if the nights
        # that this pixel is visible crosses over the 364-0 night
        # boundary

        unobservedTime = (dTs - cutoffTime).clip(0).sum()

        # These lines add unobserved time between startTime and
        # firstObs time and between lastObsTime and endTime.
        # However, the field may have been unobserved during these
        # times because they had not risen for the year (or had set)
        # and it's harder to figure this out than it is harmful to
        # just ignore this potentially unobserved time imo
        if includeEnds:
            unobservedTime += firstObsTime - startTime
            unobservedTime += max(0, endTime - lastObsTime - cutoffTime)

        return unobservedTime / (endTime - startTime - invisibleTime)

