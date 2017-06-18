from __future__ import division

import numpy as np

from lsst.sims.speedObservatory import Telescope
from lsst.sims.speedObservatory import sky
import Config
from SkyMap import SkyMap

import pygame
from matplotlib.colors import hsv_to_rgb

class GraphicalMonitor:

    def __init__(self, skyMap, mode="nvisits"):
        self.mode = mode
        supportedModes = set(["nvisits","filters"])
        if mode not in supportedModes:
            raise ValueError("graphics mode " + str(mode) + " not supported")

        # resolution is (xNumPix, yNumPix)
        resolution = skyMap.getResolution()
        # just an estimate of what will look good
        contourResolution = resolution[0] / 100

        # calculate a contour around the whole sky
        dec = np.linspace(-np.pi/2, np.pi/2, num=700*contourResolution)

        ra = np.radians(180) * np.ones(dec.shape)
        plus = skyMap.radec2imdata(np.vstack([ra, dec]).T)

        ra = np.radians(-180) * np.ones(dec.shape)
        minus = skyMap.radec2imdata(np.vstack([ra, dec]).T)

        self.wholeSkyContour = np.vstack([plus, minus])

        # calculate a contour in imdata representing 2 airmasses
        nAirmass = 2
        airmassContourAz = np.linspace(0,2 * np.pi, num=500*contourResolution)
        # airmass = csc theta
        airmassContourAlt = np.ones(airmassContourAz.shape) * np.arcsin(1 / nAirmass) 
        airmassContourAltaz = np.vstack([airmassContourAlt, airmassContourAz]).T
        self.airmassContour = skyMap.altaz2imdata(airmassContourAltaz)
        
        telescope = Telescope()

        # calculate a contour in imdata representing the zenith avoidance zone
        zenithAvoidAz = np.linspace(0, 2 * np.pi, num=30*contourResolution)
        zenithAvoidAlt = telescope.maxAlt * np.ones(zenithAvoidAz.shape)
        zenithAvoidAltaz = np.vstack([zenithAvoidAlt, zenithAvoidAz]).T
        self.zenithAvoidContour = skyMap.altaz2imdata(zenithAvoidAltaz)

        # calculate contours at various azimuths from zenith to the horizon
        azes = np.linspace(0, 2 * np.pi, num=16)
        alts = np.linspace(0, np.pi/2, num=200 * contourResolution)
        azContourAltAz = np.array([[alt, az] for az in azes for alt in alts])
        self.azContour = skyMap.altaz2imdata(azContourAltAz)

        # calculate where the meridian contour should go
        self.meridianX = skyMap.altaz2imdata(np.array([[0,0]]))[0,0]
        
        # create a lookup table for hsv=>packed color int since hsv_to_rgb is "slow"
        hue = np.linspace(0,1,num=256)
        sat = np.ones(hue.shape)
        val = np.ones(hue.shape)
        hsv = np.array([hue, sat, val]).transpose()
        rgb = (hsv_to_rgb(hsv) * 255).astype(int)
        self.packedColorLookup = (rgb[:,0] << 16) + (rgb[:,1] << 8) + rgb[:,2]

        self.filterColorLookup = np.zeros(len(Telescope.filters) + 1)
        self.filterColorLookup[-1] = (0xFF << 16) + (0xFF << 8) + 0xFF
        for i, filter in enumerate(Telescope.filters):
            if filter == "u":
                self.filterColorLookup[i] = (0x51 << 16) + (0x83 << 8) + 0xB5
            if filter == "g":
                self.filterColorLookup[i] = (0xB8 << 16) + (0x57 << 8) + 0x50
            if filter == "r":
                self.filterColorLookup[i] = (0x97 << 16) + (0xB1 << 8) + 0x5C
            if filter == "i":
                self.filterColorLookup[i] = (0x82 << 16) + (0x6A << 8) + 0x9F
            if filter == "z":
                self.filterColorLookup[i] = (0x4D << 16) + (0xA6 << 8) + 0xBC
            if filter == "y":
                self.filterColorLookup[i] = (0xE2 << 16) + (0x94 << 8) + 0x46

        self.screen = pygame.display.set_mode(resolution)

    def __del__(self):
        pygame.display.quit()

    def updateDisplay(self, skyMap, curTime):
        startTime = Config.surveyStartTime

        # skyAngle is how far to rotate the sky to the right, so we
        # need it to be -\Delta meridian since the ra of the meridian
        # increases with time (RA is higher in the East)
        skyAngle = (sky.raOfMeridian(startTime) -
                    sky.raOfMeridian(curTime))
        skyAngle = (skyAngle - skyMap.ra0) % (2*np.pi)

        if self.mode == "nvisits":
            imdata = skyMap.getNVisitsMap(skyAngle)

            if np.max(imdata) > 0:
                imdata = imdata / np.max(imdata)
            # 240/360 was chosen so we go from blue to red
            hue = (((240 - imdata * 240) / 360) * 255).astype(int).transpose()
            imdata = self.packedColorLookup[hue]

        elif self.mode == "filters":
            filterMap = skyMap.getFiltersMap(skyAngle).T
            imdata = self.filterColorLookup[filterMap]

        # draw a line down the middle to represent the meridian
        imdata[self.meridianX,:] = 0

        # outline the whole sky
        imdata[self.wholeSkyContour[:,0], self.wholeSkyContour[:,1]] = 0

        # draw a circle at zenith
        imdata[self.zenithAvoidContour[:,0], self.zenithAvoidContour[:,1]] = 0

        # draw contours of constant azimuth
        #imdata[self.azContour[:,0], self.azContour[:,1]] = 0

        # draw the contour at 2 airmasses
        imdata[self.airmassContour[:,0], self.airmassContour[:,1]] = 0

        pygame.surfarray.blit_array(self.screen, imdata)
        pygame.display.flip()
        

    def saveFrame(self, filename):
        # save the current frame to disk
	pygame.image.save(self.screen, filename)

