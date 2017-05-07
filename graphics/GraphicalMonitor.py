from __future__ import division

import numpy as np

from Telescope import Telescope
import AstronomicalSky
import Config
from minis.SkyMap import SkyMap

import pygame
from matplotlib.colors import hsv_to_rgb

class GraphicalMonitor:

    def __init__(self, skyMap):
        self.pendingVisits = []

        # calculate ra_0: the ra of the zenith at Config.surveyStartTime
        # this is used to make the zenith centered in the displayed image
        self.ra0 = AstronomicalSky.altaz2radec(np.array([[np.pi/2, 0]]),
                                               Config.surveyStartTime)[0,0]

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

        self.screen = pygame.display.set_mode(resolution)

    def updateDisplay(self, skyMap, curTime):
        startTime = Config.surveyStartTime

        # skyAngle is how far to rotate the sky to the right, so we
        # need it to be 2pi - meridian since the ra of the meridian
        # increases with time
        skyAngle = 2*np.pi - AstronomicalSky.raOfMeridian(curTime - startTime)

        nVisitsMap = skyMap.getNVisitsMap(skyAngle)
        imdata = nVisitsMap
        if np.max(nVisitsMap) > 0:
            imdata = imdata / np.max(nVisitsMap)

        # 240/360 was chosen so we go from blue to red
        hue = (((240 - imdata * 240) / 360) * 255).astype(int).transpose()
        imdata = self.packedColorLookup[hue]

        # draw a line down the middle to represent the meridian
        imdata[self.meridianX,:] = 0

        # outline the whole sky
        imdata[self.wholeSkyContour[:,0], self.wholeSkyContour[:,1]] = 0

        # draw a circle at zenith
        imdata[self.zenithAvoidContour[:,0], self.zenithAvoidContour[:,1]] = 0

        # draw contours of constant azimuth
        imdata[self.azContour[:,0], self.azContour[:,1]] = 0

        # draw the contour at 2 airmasses
        imdata[self.airmassContour[:,0], self.airmassContour[:,1]] = 0

        pygame.surfarray.blit_array(self.screen, imdata)
        pygame.display.flip()
        

    def saveFrame(self, filename):
        # save the current frame to disk
	pygame.image.save(self.screen, filename)

