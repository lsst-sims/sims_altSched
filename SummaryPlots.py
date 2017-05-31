from __future__ import division

from matplotlib import pyplot as plt
from Telescope import Telescope
import numpy as np
import AstronomicalSky as AS
import time
import itertools

import logging
logging.basicConfig(format="%(asctime)s: %(message)s")
log = logging.getLogger("SummaryPlots")
log.setLevel(logging.DEBUG)

class SummaryPlots:
    def __init__(self, skyMap, slewTimes=None):
        print
        log.debug("start __init__")

        self.skyMap = skyMap
        self.slewTimes = slewTimes

        visitInfoMap = skyMap.getVisitInfoMap()

        log.debug("got visitInfoMap")
 
        # flatten turns 2D into 1D, itertools adds all the visitInfo arrays
        # corresponding to each pixel
        allVisitInfos = visitInfoMap.flatten()
        log.debug("flattened visitInfoMap")
        allVisitInfos = list(itertools.chain.from_iterable(allVisitInfos))
        log.debug("combined all visitInfos")

        # this gives an Nx4 array (assuming visitInfo is [time, ra, dec, filter])
        # TODO make this work, and this should apply everywhere
        # dtype = [("time",float), ("ra", float), ("dec", float), ("filter", "S1")]
        #allVisitInfos = np.array(map(tuple, allVisitInfos), dtype=np.dtype(dtype))
        allVisitInfos = np.array(allVisitInfos)

        log.debug("calculated allVisitInfos")

        times = allVisitInfos[:,0].astype(float)
        self.ras = allVisitInfos[:,1].astype(float)
        self.decs = allVisitInfos[:,2].astype(float)
        self.filters = allVisitInfos[:,3].astype("S1")

        log.debug("sliced allVisitInfos")
        altazes = AS.radec2altaz(np.vstack([self.ras, self.decs]).T, times)
        self.alts = altazes[:,0]
        self.azes = altazes[:,1]#.clip(0,2*np.pi)

        log.debug("got altazes")

        # now get all revisit times
        # (one for loop is necessary because the number of visitInfos in a
        #  pixel is variable)
        def getRevisitTimes(pix):
            revisitTimes = []
            prevVisitTime = -3600*24*365
            for visitInfo in pix:
                time = visitInfo[0]
                revisitTime = time - prevVisitTime
                if revisitTime < 3600 * 24 * 30 * 4:
                    revisitTimes.append(revisitTime)
                prevVisitTime = time
            return revisitTimes
        self.revisitTimes = np.vectorize(getRevisitTimes, otypes="O")(visitInfoMap)
        self.revisitTimes = self.revisitTimes.flatten()

        log.debug("got revisitTimes")

        self.revisitTimes = list(itertools.chain.from_iterable(self.revisitTimes))

        log.debug("flattened revisitTimes")
        log.debug("end __init__")

    def show(self):
        # called to show plots that have been queued
        plt.show()
 
    def dAirmassCum(self):
        # calculate dAirmass for each visit
        dAirmasses = []
        log.debug("start dAirmassCum")

        for alt, dec in zip(self.alts, self.decs):
            airmass = 1 / np.sin(alt)
            dAirmass = airmass - (1/np.cos(dec - Telescope.latitude))
            dAirmasses.append(dAirmass)

        log.debug("got dAirmasses")

        # make a histogram showing the cumulative distribution of dAirmasses
        plt.figure("DAirmass Cumulative Distribution")
        plt.title("Cumulative distribution of observed airmass - optimal airmass")
        plt.xlabel("Number of airmasses observed through - " + \
                       "sec(declination - telescope latitude)")
        plt.ylim(0,1)
        plt.hist(dAirmasses, 300, range=(0,0.6), cumulative=True, normed=True)
        log.debug("end dAirmassCum")
    
    def airmassHist(self):
        airmasses = []
        for alt in self.alts:
            airmasses.append(1 / np.sin(alt))
        # now make a histogram of the airmass observed through
        plt.figure("Airmass Histogram")
        plt.title("Histogram of airmasses observed through")
        plt.xlabel("airmass = csc(altitude)")
        plt.hist(airmasses, 300)

    def dAirmassContour(self):
        # dAirmass = delta airmass = observed airmass - optimal airmass given dec
        # first, create a polar plot with r = dAirmass and theta = azimuth

        # limits and resolution for dAirmass
        nDAirmass = 50
        minDAirmass = 0
        maxDAirmass = 0.5

        # limits for az are always 0 to 2pi
        # doesn't make sense to have resolution higher than ~3 degrees
        # since the FoV size will smear out the graph on that scale
        nAzes = 100

        # r, theta are polar coordinates in the plot
        r, theta = np.meshgrid(np.linspace(minDAirmass, maxDAirmass, num=nDAirmass),
                               np.linspace(0, 2*np.pi, num=nAzes))

        values = np.zeros(r.shape)

        for alt, az, dec in zip(self.alts, self.azes, self.decs):
            # airmass number is csc of elevation or sec of zenith angle
            airmass = 1 / np.sin(alt)

            # the best place to observe a field is when it hits the meridian,
            # at which point you are observing through an airmass of
            # sec(dec - latitude)
            dAirmass = airmass - (1/np.cos(dec - Telescope.latitude))
            
            # don't plot points out of range
            if dAirmass > maxDAirmass or dAirmass < minDAirmass:
                continue

            # rId and thetaId are indices within r and theta corresponding
            # to the dAirmass/zenith and azimuth values
            rId = (dAirmass - minDAirmass) / (maxDAirmass - minDAirmass) * nDAirmass
            rId = int(rId)

            thetaId  = int(az / (2*np.pi) * (nAzes-1))
            values[thetaId, rId] += 1

        # now show contour plots in polar coordinates where r is dAirmass
        # or zenith angle and theta is azimuth
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.set_title("log_10(nVisits)")
        ax.set_theta_zero_location("N")
        ax.set_ylim(minDAirmass, maxDAirmass)
        ax.set_xlabel("Azimuth (N=0deg, E=90deg)")
        ax.set_ylabel("Delta Airmass (obsAirmass - sec(obsDec - Tel.lat))")
        p = ax.contourf(theta, r, np.log10(values))
        plt.colorbar(p, ax=ax)

    def zenithAngleContour(self):
        # similar to self.dAirmassContour
        log.debug("start zenithAngleContour")

        nZenith = 75
        minZenith = 0
        maxZenith = 90

        nAzes = 150

        # create a grid of r and theta
        r, theta = np.meshgrid(np.linspace(minZenith, maxZenith, num=nZenith),
                               np.linspace(0, 2*np.pi, num=nAzes))
        values = np.zeros(r.shape)
        for alt, az, dec in zip(self.alts, self.azes, self.decs):
            zenithAngle = 90 - np.degrees(alt)
            if zenithAngle < minZenith or zenithAngle >= maxZenith:
                continue

            # calculate offsets into r and theta
            rId = (zenithAngle - minZenith) / (maxZenith - minZenith) * nZenith
            rId = int(rId)
            thetaId = int(az / (2*np.pi) * (nAzes-1))
            values[thetaId, rId] += 1

        log.debug("got values")

        # and plot
        fix, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.set_title("nVisits")
        ax.set_theta_zero_location("N")
        ax.set_ylim(minZenith, maxZenith)
        ax.set_xlabel("Azimuth (N=0deg, E=90deg)")
        ax.set_ylabel("Zenith angle (degrees)")
        p = ax.contourf(theta, r, values)
        plt.colorbar(p, ax=ax)

        log.debug("end zenithAngleContour")

    def slewHist(self, maxSlew = 30):
        plt.figure("Slew Time Histogram")

        histBins = np.arange(min(self.slewTimes), maxSlew, 0.1)
        plt.hist(self.slewTimes, bins = histBins)
        plt.xlabel("Slew Time (secs)")
        plt.ylabel("Number of slews")
        plt.title("Histogram of Slew Times")

    def slewRankCum(self):
        plt.figure("Slew Time Cumulative Rank Distribution")

        plt.title("Cumulative sum of slew times")
        plt.xlabel("Cumulative sum up to the nth fastest slew")
        plt.ylabel("Cumulative sum (secs)")
        sortedTimes = np.sort(self.slewTimes)
        cum = np.cumsum(sortedTimes)
        plt.plot(np.arange(len(cum)), cum)

    def revisitHist(self):
        revisitTimesMins = np.array(self.revisitTimes) / 60
        plt.figure("Revisit Time Histogram")
        plt.hist(revisitTimesMins, 
                 bins = np.arange(0, 3*np.median(revisitTimesMins), 1))
        plt.xlabel("Revisit Time (mins)")
        plt.ylabel("Number of Revisits")
        plt.title("Histogram of Revisit Times")
