from __future__ import division

from astropy import wcs
import numpy as np

import AstronomicalSky
import Telescope
import Utils
import Config

import palpy

class SkyMap:
    def radec2imdata(self, radec):
        """
        This method takes in pixels in the form (RA, DEC) and returns
        indices into the imdata array (y, x) 
        """
        # * subtract off ra_0 so ra=ra_0 is centered instead of ra=0
        # * convert from radec to raw wcs pixels
        # * multiply by resScale
        # * multiply dec by -1 so higher dec => higher in the img, not lower
        # * center the result 
        
        rawPix = self.w.wcs_world2pix(np.degrees(radec),1)
        imdataPix = rawPix * self.resScale
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
        # * divide by resScale
        # * convert the resulting raw wcs pixels to radec
        
        rawPix = imdata - np.array([self.xCenter, self.yCenter])
        rawPix[:,0] *= -1
        rawPix /= self.resScale
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
        pix *= self.resScale

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

    def __init__(self, resScale):
        self.resScale = resScale 
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
        (self.yMin, self.yMax) = (int(resScale * wcsYMin), int(resScale * wcsYMax))
        (self.xMin, self.xMax) = (int(resScale * wcsXMin), int(resScale * wcsXMax))
        self.xCenter = int((self.xMax - self.xMin) / 2)
        self.yCenter = int((self.yMax - self.yMin) / 2)

        self.ra0 = AstronomicalSky.altaz2radec(np.array([[np.pi/2, 0]]),
                                               Config.surveyStartTime)[0,0]
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
        self.skyPix = self.w.wcs_pix2world(projPixReversed / resScale, 1)
        self.skyPix = np.radians(self.skyPix)
        print "skyPix", self.skyPix


        # take out pixels that are in the projection rectangle but
        # don't actually represent places on the sky (i.e. pixels in the corners)
        validMask = ~np.any(np.isnan(self.skyPix), axis=1)
        projPix = projPix[validMask]
        print "projPix", projPix.shape
        self.skyPix = self.skyPix[validMask]
        print "self.skyPix[validMask]", self.skyPix, self.skyPix.shape

        # store the number of visits to each pixel in pixValues
        self.pixValues = np.zeros(len(self.skyPix))

        # keep track of how long each row of pixels is in the mollweide
        self.rowLengths = {y: (projPix[:,0] == y).sum() 
                           for y in range(self.yMin, self.yMax)}
        self.rowLengthsCumSum = [sum(self.rowLengths[y] 
                                   for y in range(self.yMin, self.yMin + i))
                                   for i in range(self.yMax - self.yMin)]
        self.rowLengthsCumSum = np.array(self.rowLengthsCumSum)

    def getResolution(self):
        return (self.xMax - self.xMin, self.yMax - self.yMin)
 
    def clear(self):
        # clear pixValues
        self.pixValues = np.zeros(self.pixValues.shape)

    def addVisit(self, visit, time):
        # edgeOfFov should be the radius of a circle that includes all of the fov
        edgeOfFov = Utils.spherical2Cartesian(0, 1.5 * Telescope.fovWidth / 2)
        r = np.linalg.norm(np.array([1,0,0]) - np.array(edgeOfFov))

        # get a list of all pixels which might lie on the focal plane
        # (I was originally using a KDTree of all sky pixels in cartesian
        # coordinates and then querying it to get the candidate pix, but
        # this was about 20x slower (for resScale=5) )

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
        candidatePixIds = candidatePixIds[withinRange]

        candidateSkyPix = self.skyPix[candidatePixIds]

        # calculate the pupil coordinates x, y of the obsMetaData
        # i.e. project candidateSkyPix onto the plane tangent to the unit
        # sphere at the altaz of (visit.ra, visit.dec)

        radec = np.array([[visit.ra, visit.dec]])
        visitAltaz = AstronomicalSky.radec2altaz(radec, time)
        candidateAltaz = AstronomicalSky.radec2altaz(candidateSkyPix, time)

        x, y = palpy.ds2tpVector(candidateAltaz[:,1], candidateAltaz[:,0],
                                 visitAltaz[0][1], visitAltaz[0][0])
        x *= -1

        # TODO rotate depending on the angle of the focal plane wrt the sky
        # I think I need to add/subtract the parallactic angle?
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
        self.pixValues[coveredPixIds] += 1

    def get2DMap(self, skyAngle):
        # this array will hold the 2D, rotated image
        imdata = np.zeros((self.yMax - self.yMin, self.xMax - self.xMin))

        # the mollweide projection has increasing RA to the right
        # which means that we're looking at the front of the celestial
        # sphere with earth behind. The celestial sphere rotates clockwise 
        # so the projection should move to the left as time increases
        # I believe this goes against convention but too bad

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
            # TODO I got an error here once:
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

        #print np.max(self.pixValues)
        imdata /= np.max(self.pixValues)

        return imdata

