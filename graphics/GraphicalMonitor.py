from __future__ import division

import numpy as np
from numpy import cross, eye, dot
from scipy.linalg import expm3, norm
from scipy.spatial import KDTree
from astropy import wcs

import matplotlib
from matplotlib import pyplot as plt
from multiprocessing import Pool
from collections import deque
from collections import Counter
import itertools

import Telescope
import Utils
import AstronomicalSky
import Config


import palpy

import pygame
from matplotlib.colors import hsv_to_rgb


class GraphicalMonitor:

    def radec2imdata(self, radec):
        """
        This method takes in pixels in the form (RA, DEC) and returns
        indices into the imdata array (y, x) 
        """
        # * subtract off ra_0 so ra=ra_0 is centered instead of ra=0
        # * convert from radec to raw wcs pixels
        # * multiply by imScale
        # * multiply dec by -1 so higher dec => higher in the img, not lower
        # * center the result 
        
        rawPix = self.w.wcs_world2pix(np.degrees(radec),1)
        imdataPix = rawPix * self.imScale
        imdataPix[:,1] *= -1
        imdataPix = imdataPix.astype(int)
        imdataPix += np.array([self.xCenter,self.yCenter])
        return imdataPix

    def imdata2radec(self, imdata):
        """
        Converts indices into the imdata array into radec 
        """
        # * subtract out the center
        # * multiply y by -1
        # * divide by imScale
        # * convert the resulting raw wcs pixels to radec
        
        rawPix = imdata - np.array([self.xCenter, self.yCenter])
        rawPix[:,0] *= -1
        rawPix /= self.imScale
        radec = self.w.wcs_pix2world(np.degrees(rawPix),1)
        return radec

    def altaz2imdata(self, altaz):
        """
        Converts alt/az to indices into the imdata array at the time
        specified by Config.surveyStartTime
        """
        radec = AstronomicalSky.altaz2radec(altaz, Config.surveyStartTime)
        radec[:,0] -= self.ra0
        return self.radec2imdata(radec)

    def imdata2altaz(self, imdata):
        """
        Converts indices into the imdata array into alt/az at the time
        specified by Config.surveyStartTime
        """
        raise NotImplementedError("imdata2altaz not implemented yet")
        return

    def radec2projPix(self, radec):
        # pix is the coord where the center is (0,0)
        pix = self.w.wcs_world2pix(np.degrees(radec), 1)
        pix *= self.imScale

        # ys is the y value where 0 is in the middle
        ys = pix[:,1].astype(int)

        # rowNums is the y value where 0 is at the top
        rowNums = ys - self.yMin
        rowLengths = np.array([self.rowLengths[y] for y in ys])

        # rowStarts is the index into self.skyPix of the start
        # of the row that radec is in
        rowStarts = self.rowLengthsCumSum[rowNums]

        # offsets is the x offset needed to add to rowStarts
        # to get the actual position of pix within self.skyPix
        offsets = (pix[:,0] + (rowLengths / 2)).astype(int)

        return rowStarts + offsets

    def __init__(self, context, imScale, numVisitsDisplayed=None):
        self.context = context
        self.imScale = imScale
        self.pendingVisits = []

        # create a WCS object to do the transforms from the sphere to the
        # Mollweide projection
        self.w = wcs.WCS(naxis=2)
        # TODO figure out what the dash syntax is
        self.w.wcs.ctype = ["RA---MOL", "DEC--MOL"]

        # this the range of output from wcs_world2pix
        yLim = np.ceil(self.w.wcs_world2pix([[0,90]], 1)[0,1])
        xLim = np.ceil(self.w.wcs_world2pix([[180,0]], 1)[0,0])
        (wcsYMin, wcsYMax) = (-yLim, yLim)
        (wcsXMin, wcsXMax) = (-xLim, xLim)

        # these are the range of pixels we're going to display
        (self.yMin, self.yMax) = (int(imScale * wcsYMin), int(imScale * wcsYMax))
        (self.xMin, self.xMax) = (int(imScale * wcsXMin), int(imScale * wcsXMax))
        self.xCenter = int((self.xMax - self.xMin) / 2)
        self.yCenter = int((self.yMax - self.yMin) / 2)

        # a list of all the pixels [y, x]
        # note that the order that these are created in is very important
        # since the innter loop in updateDisplay assumes this ordering 
        projPix = np.array([[y, x] for y in range(self.yMin, self.yMax)
                            for x in range(self.xMin, self.xMax)])

        # convert the projected pixels to sky coordinates
        # TODO I think maybe the origin should be 0 not 1? 
        # The ra values only go up to 359
        # TODO refactor
        projPixReversed = np.vstack([projPix[:,1], projPix[:,0]]).T
        self.skyPix = self.w.wcs_pix2world(projPixReversed / imScale, 1)
        self.skyPix = np.radians(self.skyPix)
        print "skyPix", self.skyPix

        # calculate ra_0: the ra of the zenith at Config.surveyStartTime
        # this is used to make the zenith centered in the displayed image
        self.ra0 = AstronomicalSky.altaz2radec(np.array([[np.pi/2, 0]]),
                                               Config.surveyStartTime)[0,0]

        # calculate a contour around the whole sky
        dec = np.linspace(-np.pi/2, np.pi/2, num=700*self.imScale)

        ra = np.radians(180) * np.ones(dec.shape)
        plus = self.radec2imdata(np.vstack([ra, dec]).T)

        ra = np.radians(-180) * np.ones(dec.shape)
        minus = self.radec2imdata(np.vstack([ra, dec]).T)

        self.wholeSkyContour = np.vstack([plus, minus])

        # calculate a contour in imdata representing 2 airmasses
        nAirmass = 2
        airmassContourAz = np.linspace(0,2 * np.pi, num=500*self.imScale)
        # airmass = csc theta
        airmassContourAlt = np.ones(airmassContourAz.shape) * np.arcsin(1 / nAirmass) 
        airmassContourAltaz = np.vstack([airmassContourAlt, airmassContourAz]).T
        self.airmassContour = self.altaz2imdata(airmassContourAltaz)
        
        # calculate a contour in imdata representing the zenith avoidance zone
        zenithAvoidAz = np.linspace(0, 2 * np.pi, num=30*self.imScale)
        zenithAvoidAlt = Telescope.maxAlt * np.ones(zenithAvoidAz.shape)
        zenithAvoidAltaz = np.vstack([zenithAvoidAlt, zenithAvoidAz]).T
        self.zenithAvoidContour = self.altaz2imdata(zenithAvoidAltaz)

        # calculate contours at various azimuths from zenith to the horizon
        azes = np.linspace(0, 2 * np.pi, num=16)
        alts = np.linspace(0, np.pi/2, num=200 * self.imScale) 
        azContourAltAz = np.array([[alt, az] for az in azes for alt in alts])
        self.azContour = self.altaz2imdata(azContourAltAz)

        # calculate where the meridian contour should go
        self.meridianX = self.altaz2imdata(np.array([[0,0]]))[0,0]

        # take out pixels that are in the projection rectangle but
        # don't actually represent places on the sky (i.e. pixels in the corners)
        validMask = ~np.any(np.isnan(self.skyPix), axis=1)
        projPix = projPix[validMask]
        print "projPix", projPix.shape
        self.skyPix = self.skyPix[validMask]
        print "self.skyPix[validMask]", self.skyPix, self.skyPix.shape

        # put the sky coordinates into a k-d tree so we can quickly 
        # compute nearest neighbors
        cartesianPix = Utils.spherical2Cartesian(self.skyPix[:,0], self.skyPix[:,1])
        print "cartesianPix", cartesianPix
        self.skyPixTree = KDTree(cartesianPix)

        # store the number of visits to each pixel in pixValues
        self.pixValues = np.zeros(len(self.skyPix))
        
        if numVisitsDisplayed is None:
            self.useVisitQueue = False
        else:
            self.useVisitQueue = True
            self.visitQueue = deque(maxlen=numVisitsDisplayed)

        # keep track of how long each row of pixels is in the mollweide
        self.rowLengths = {y: (projPix[:,0] == y).sum() 
                           for y in range(self.yMin, self.yMax)}
        self.rowLengthsCumSum = [sum(self.rowLengths[y] 
                                   for y in range(self.yMin, self.yMin + i))
                                 for i in range(self.yMax - self.yMin)]
        self.rowLengthsCumSum = np.array(self.rowLengthsCumSum)
        
        # create a lookup table for hsv=>packed color int since hsv_to_rgb is "slow"
        hue = np.linspace(0,1,num=256)
        sat = np.ones(hue.shape)
        val = np.ones(hue.shape)
        hsv = np.array([hue, sat, val]).transpose()
        rgb = (hsv_to_rgb(hsv) * 255).astype(int)
        self.packedColorLookup = (rgb[:,0] << 16) + (rgb[:,1] << 8) + rgb[:,2]

        self.screen = pygame.display.set_mode((self.xMax - self.xMin, 
                                               self.yMax - self.yMin))


    def addVisit(self, visit):
        self.pendingVisits.append(visit)

    def _addVisits(self):
        # batch process the pending visits added using self.addVisit
        # edgeOfFov should be the radius of a circle that includes all of the fov
        edgeOfFov = Utils.spherical2Cartesian(0, 1.5 * Telescope.fovWidth / 2)
        r = np.linalg.norm(np.array([1,0,0]) - np.array(edgeOfFov))

        ras =  np.array([visit.ra  for visit in self.pendingVisits])
        decs = np.array([visit.dec for visit in self.pendingVisits])

        # N x 3
        cartesianVisits = Utils.spherical2Cartesian(ras, decs)

        # get a list of all pixels which might lie on the focal plane
        # (I was originally using a KDTree of all sky pixels in cartesian
        # coordinates and then querying it to get the candidate pix, but
        # this was about 20x slower (for imScale=5) )

        for i, visit in enumerate(self.pendingVisits):
            # candidate pix must be within 2 degrees of visit dec
            deltaDec = np.radians(2)
            minDec = max(visit.dec - deltaDec, -np.pi/2)
            maxDec = min(visit.dec + deltaDec, np.pi/2)
            projPixRange = self.radec2projPix([[0, minDec],
                                               [0, maxDec]])
            
            candidatePixIds = np.arange(projPixRange[0], projPixRange[1])

            candidateRas, candidateDecs = self.skyPix[candidatePixIds].T

            # candidate pix must be within 2 degrees of visit ra normalized
            # by cos(dec) to avoid cutting off near poles
            # np.max to prevent div by 0
            deltaRa = np.radians(2) / max(np.cos(candidateDecs).max(), 0.001)
            
            raRange = ((visit.ra - deltaRa) % (2*np.pi), 
                       (visit.ra + deltaRa) % (2*np.pi))
            withinRange = Utils.areRasInRange(candidateRas, raRange)
            #withinRange = np.where(abs(candidateRas - visit.ra) < deltaRa)
            candidatePixIds = candidatePixIds[withinRange]

            candidateSkyPix = self.skyPix[candidatePixIds]

            # calculate the pupil coordinates x, y of the obsMetaData
            # i.e. project candidateSkyPix onto the plane tangent to the unit
            # sphere at the altaz of (visit.ra, visit.dec)

            radec = np.array([[visit.ra, visit.dec]])
            visitAltaz = AstronomicalSky.radec2altaz(radec, self.context.time())
            candidateAltaz = AstronomicalSky.radec2altaz(candidateSkyPix,
                                                         self.context.time())

            x, y = palpy.ds2tpVector(candidateAltaz[:,1], candidateAltaz[:,0],
                                     visitAltaz[0][1], visitAltaz[0][0])
            x *= -1

            # TODO rotate depending on the angle of the focal plane wrt the sky
            # not just the absolute rotation of the visit
            xPupil = x * np.cos(visit.rotation) - y * np.sin(visit.rotation)
            yPupil = x * np.sin(visit.rotation) + y * np.cos(visit.rotation)

            # figure out which chip these pupil coords lie on
            pupilCoords = np.vstack([xPupil, yPupil]).T

            def arePixCovered(x, y):
                # this code ignores the gaps between rafts and chips

                # I found these by evaluating chilNameFromPupilCoords many times
                # please ignore the number of sig figs (no idea if my testing
                # was that accurate)
                # TODO I think this might be in units of radians?
                edge = 0.030756 # outer edge in units of pupil coords
                raftWidth = 0.012303 # width of a raft in pupil coords

                # consider the camera CCD collection to be a cross composed
                # of a horizontal and vertical rectangle placed on top of
                # each other. The first two lines check if the pixel is in
                # the horizontal rectangle and the second two in the vertical

                return (((x > -edge) & (x < edge) &
                         (y > -edge + raftWidth) & (y < edge - raftWidth)) |
                        ((y > -edge) & (y < edge) &
                         (x > -edge + raftWidth) & (x < edge - raftWidth)))

            
            coveredPixIds = candidatePixIds[arePixCovered(pupilCoords[:,0],
                                                          pupilCoords[:,1])]
            if self.useVisitQueue:
                self.visitQueue.append(coveredPixIds)
            else:
                # if we're not keeping the visits in a queue, we can just 
                # increment coveredPixIds's counts in pixValues
                self.pixValues[coveredPixIds] += 1

        # reset pendingVisits
        self.pendingVisits = []


    def updateDisplay(self):
        self._addVisits()

        startTime = Config.surveyStartTime
        curTime = self.context.time()

        if self.useVisitQueue:
            # construct pixValues from self.visitQueue
            visitedPix = list(itertools.chain.from_iterable(self.visitQueue))
            pixCounts = Counter(visitedPix)
            self.pixValues = np.zeros(self.pixValues.shape)
            self.pixValues[pixCounts.keys()] = pixCounts.values()
        else:
            # the correct counts are already in self.pixValues
            pass
        

        # this array will hold the rotated image
        imdata = np.zeros((self.yMax - self.yMin, self.xMax - self.xMin))

        # the mollweide projection has increasing RA to the right
        # which means that we're looking at the front of the celestial
        # sphere with earth behind. The celestial sphere rotates clockwise 
        # so the projection should move to the left as time increases
        # I believe this goes against convention but too bad

        # skyAngle is how far to rotate the sky to the right, so we
        # need it to be 2pi - meridian since the ra of the meridian
        # increases with time
        skyAngle = 2*np.pi - AstronomicalSky.raOfMeridian(curTime - startTime)

        pixMvmtRatio = skyAngle / (2*np.pi)
        
        # keep track of where in pixValues the start of the current row is
        rowStartProjPixId = 0
        for y in xrange(self.yMin, self.yMax):
            # y is the y coordinate in the image (so yMax is the top row)

            # iy is the y index into imdata (so 0 is the top row)
            iy = (self.yMax - 1) - y

            rowLen = self.rowLengths[y]

            # these are indices into imdata[y], not into pixvalues
            rowStart = int((self.xMax - self.xMin - rowLen) / 2)
            rowEnd = rowStart + rowLen

            # calculate how many pixels to shift this row by
            numPixMoved = int(rowLen * pixMvmtRatio)

            
            # first handle the pixels in this row which will be replaced by
            # pixels which have wrapped around
            # TODO I get an error here:
            # "could not broadcast input array from shape (0) into shape (662)"
            # around day 125
            imdata[iy, rowStart:rowStart + numPixMoved] = \
                    self.pixValues[rowStartProjPixId + rowLen - numPixMoved : \
                                   rowStartProjPixId + rowLen]

            # then handle the rest of the pixels 
            imdata[iy, rowStart + numPixMoved:rowStart + rowLen] = \
                    self.pixValues[rowStartProjPixId: \
                                   rowStartProjPixId + rowLen - numPixMoved]

            # keep track of where in pixValues the next row starts
            rowStartProjPixId += rowLen

        imdata /= np.max(self.pixValues)

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


