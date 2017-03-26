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

#from lsst.sims.coordUtils import _chipNameFromRaDec
#from lsst.sims.coordUtils import chipNameFromPupilCoords
#from lsst.sims.utils import ObservationMetaData
#from lsst.obs.lsstSim import LsstSimMapper

import palpy

import pygame
from matplotlib.colors import hsv_to_rgb

# get the LSST camera so we can figure out which pixels are covered
# by the camera's sensors
#mapper = LsstSimMapper()
#camera = mapper.camera

def workerWrapper(pupilCoords):
    return chipNameFromPupilCoords(pupilCoords[:,0], 
                                   pupilCoords[:,1], 
                                   camera = camera)

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
        
        #shifted = radec
        #shifted[:,0] -= self.ra0
        shifted = radec
        rawPix = self.w.wcs_world2pix(np.degrees(shifted),0)
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
        radec = self.w.wcs_pix2world(np.degrees(rawPix),0)
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

    def __init__(self, context, imScale, numVisitsDisplayed=None):
        self.context = context
        self.imScale = imScale

        # this the approximate range of output from wcs_world2pix
        (wcsYMin, wcsYMax) = (-85, 85)
        (wcsXMin, wcsXMax) = (-170, 170)

        # these are the range of pixels we're going to display
        (self.yMin, self.yMax) = (int(imScale * wcsYMin), int(imScale * wcsYMax))
        (self.xMin, self.xMax) = (int(imScale * wcsXMin), int(imScale * wcsXMax))
        self.xCenter = int((self.xMax - self.xMin) / 2)
        self.yCenter = int((self.yMax - self.yMin) / 2)

        # create a WCS object to do the transforms from the sphere to the
        # Mollweide projection
        self.w = wcs.WCS(naxis=2)
        # TODO figure out what the dash syntax is
        self.w.wcs.ctype = ["RA---MOL", "DEC--MOL"]

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
        self.skyPix = self.w.wcs_pix2world(projPixReversed / imScale, 0)
        self.skyPix = np.radians(self.skyPix)

        """
        # figure out where the zenith is
        zenithDec = Telescope.latitude 
        zenithImdata = self.radec2imdata(np.array([[0, zenithDec]]))[0]
        zenithWorld = np.degrees(np.array([[0,zenithDec]]))
        zenithPix = w.wcs_world2pix(zenithWorld,1)
        zenithPix = self.worldPix2Imdata(zenithPix)
        """

        # calculate ra_0: the ra of the zenith at Config.surveyStartTime
        # this is used to make the zenith centered in the displayed image
        self.ra0 = AstronomicalSky.altaz2radec(np.array([[np.pi/2, 0]]),
                                               Config.surveyStartTime)[0,0]

        # calculate a contour around the whole sky
        dec = np.linspace(-np.pi/2, np.pi/2, num=2000)

        ra = np.radians(180) * np.ones(dec.shape)
        plus = self.radec2imdata(np.vstack([ra, dec]).T)

        ra = np.radians(-180) * np.ones(dec.shape)
        minus = self.radec2imdata(np.vstack([ra, dec]).T)

        self.wholeSkyContour = np.vstack([plus, minus])

        # calculate a contour in imdata representing 2 airmasses
        nAirmass = 2
        airmassContourAz = np.linspace(0,2 * np.pi, num=500)
        # airmass = csc theta
        airmassContourAlt = np.ones(airmassContourAz.shape) * np.arcsin(1 / nAirmass) 
        airmassContourAltaz = np.vstack([airmassContourAlt, airmassContourAz]).T
        self.airmassContour = self.altaz2imdata(airmassContourAltaz)
        """
        twoAirmassRaDec = AstronomicalSky.altaz2radec(twoAirmassAltaz
                                                      Config.surveyStartTime - 6.6 * 3600)
        #print "radec", twoAirmassRaDec
        twoAirmassWorld = np.vstack([twoAirmassRaDec[0], twoAirmassRaDec[1]]).T
        twoAirmassPix = w.wcs_world2pix(np.degrees(twoAirmassWorld),0)
        self.twoAirmassContour = self.worldPix2Imdata(twoAirmassPix)
        #print "twoAirmassContour", self.twoAirmassContour
        """

        """
        zenithCircleRadius = 5
        alpha = np.linspace(0,2 * np.pi, num=50)
        unitCircle = np.vstack([np.cos(alpha), np.sin(alpha)])
        self.zenithCircle = (zenithCircleRadius * unitCircle.T).astype(int)
        self.zenithCircle += zenithPix
        """

        # calculate a contour in imdata reprezenting the zenith avoidance zone
        zenithAvoidAz = np.linspace(0, 2 * np.pi, num = 100)
        zenithAvoidAlt = Telescope.maxAlt * np.ones(zenithAvoidAz.shape)
        zenithAvoidAltaz = np.vstack([zenithAvoidAlt, zenithAvoidAz]).T
        self.zenithAvoidContour = self.altaz2imdata(zenithAvoidAltaz)


        # take out pixels that are in the projection rectangle but
        # don't actually represent places on the sky (i.e. pixels in the corners)
        validMask = ~np.any(np.isnan(self.skyPix), axis=1)
        projPix = projPix[validMask]
        self.skyPix = self.skyPix[validMask]

        # put the sky coordinates into a k-d tree so we can quickly 
        # compute nearest neighbors
        cartesianPix = Utils.spherical2Cartesian(self.skyPix[:,0], self.skyPix[:,1])
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

        #self.p = Pool(8)
        
        """ for matplotlib
        # turn on interactive mode
        plt.ion()
        fig, ax = plt.subplots(figsize=(24,12))

        #ax.set_xlim(0, xMax - xMin)
        #ax.set_ylim(0, yMax - yMin)

        # show zeros to start with (we'll dynamically update it later)
        fillerImData = np.zeros((self.yMax - self.yMin, self.xMax - self.xMin))
        cmap = matplotlib.cm.jet
        cmap.set_bad("black",alpha=1)
        self.image = ax.imshow(fillerImData, vmin=0, cmap=cmap)
        plt.colorbar(self.image, shrink=0.6)
        plt.tight_layout()

        # fig.canvas.draw() works for most backends but isn't supported
        # in the one I'm using, so I need to use plt.pause() instead
        plt.pause(0.001)
        #fig.canvas.draw()
        """
        
        # create a lookup table for hsv=>packed color int since hsv_to_rgb is slow
        hue = np.linspace(0,1,num=256)
        sat = np.ones(hue.shape)
        val = np.ones(hue.shape)
        hsv = np.array([hue, sat, val]).transpose()
        rgb = (hsv_to_rgb(hsv) * 255).astype(int)
        self.packedColorLookup = (rgb[:,0] << 16) + (rgb[:,1] << 8) + rgb[:,2]

        self.screen = pygame.display.set_mode((self.xMax - self.xMin, 
                                               self.yMax - self.yMin))


    def __del__(self):
        #self.p.close()
            #self.p.join()
        pass

    def addVisit(self, visit):
        #obsMetaData = ObservationMetaData(pointingRA = np.degrees(visit.ra),
        #                                  pointingDec = np.degrees(visit.dec),
        #                                  rotSkyPos = np.degrees(visit.rotation),
        #                                  mjd = Utils.mjd(self.context.time()))
        # edgeOfFov should be the radius of a circle that includes all of the fov
        edgeOfFov = Utils.spherical2Cartesian(0, 1.5 * Telescope.fovWidth / 2)
        r = np.linalg.norm(np.array([1,0,0]) - np.array(edgeOfFov))
        cartesianVisit = Utils.spherical2Cartesian(visit.ra, visit.dec)

        # get a list of all pixels which might lie on the focal plane
        # TODO this takes ~2-10ms. Just running .query() takes 0.5ms
        # conclusion: precompute nearest neighbor balls around every pixel
        # and then here all we have to do is .query() to get the nearest
        # pixel followed by an index into the precomputed array (~0.5us)
        # looks like precomputing takes 23 seconds (I think this is for
        # imscale=2)

        candidatePixIds = self.skyPixTree.query_ball_point(cartesianVisit, r)
        candidatePixIds = np.array(candidatePixIds)
        """
        Ideally I could just use this, but _chipNameFromRaDec calls 
        _pupilCoordsFromRaDec which first converts from ICRS ra/dec to 
        observed ra/dec which is all well and good except that it's 
        painfully slow (6ms/pixel). So instead I just do what 
        _pupilCoordsFromRaDec does but cut the conversion

        chipNames = _chipNameFromRaDec(skyPix[candidatePixIds,0],
                                       skyPix[candidatePixIds,1],
                                       camera = camera,
                                       epoch = 2000.0,
                                       obs_metadata = obsMetaData)
        """
        candidateSkyPix = self.skyPix[candidatePixIds]

        # calculate the pupil coordinates x, y of the obsMetaData
        # i.e. project candidateSkyPix onto the plane tangent to the unit
        # sphere at point (visit.ra, visit.dec)
        # TODO you could probably get away with not even calling this
        # method for points near the equator given the whole small angle
        # approximation thing
        x, y = palpy.ds2tpVector(candidateSkyPix[:,0], candidateSkyPix[:,1],
                                 visit.ra, visit.dec)
        x *= -1

        # TODO rotate depending on the angle of the focal plane wrt the sky
        # not just the absolute rotation of the visit
        xPupil = x * np.cos(visit.rotation) - y * np.sin(visit.rotation)
        yPupil = x * np.sin(visit.rotation) + y * np.cos(visit.rotation)

        # figure out which chip these pupil coords lie on
        pupilCoords = np.vstack([xPupil, yPupil]).T

        def arePixCovered(pupilCoords):
            """ This is veeeery slow (~3/4 of the runtime of the entire program)
            chipNames = self.p.map(workerWrapper, np.array_split(pupilCoords, 8)) 
            chipNames = np.hstack(chipNames)
            chipNames = chipNameFromPupilCoords(xPupil, 
                                                yPupil, 
                                                camera = camera)

            # chipNameFromPupilCoords returns None if the pixel does not fall
            # on one of the camera sensors
            return np.where(~np.equal(chipNames, None))[0]
            """

            # instead, this pretty close to what the above code would return
            # but ignores the gaps between rafts and chips

            # I found these by evaluating chilNameFromPupilCoords many times
            # please ignore the number of sig figs (no idea if my testing
            # was that accurate)
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

        
        coveredPixIds = candidatePixIds[arePixCovered(pupilCoords)]
        if self.useVisitQueue:
            self.visitQueue.append(coveredPixIds)
        else:
            # if we're not keeping the visits in a queue, we can just 
            # increment coveredPixIds's counts in pixValues
            self.pixValues[coveredPixIds] += 1


    def updateDisplay(self):
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

        # values chosen so we go from blue to red
        hue = (240 - imdata * 240) / 360

        """ this does the hsv conversion manually instead of using the lookup
        # we want hsv to have shape (x, y, 3)
        hsv = np.array([hue, np.ones(hue.shape), np.ones(hue.shape)]).transpose([2,1,0])

        # hsv_to_rgb is slow, so initialize rgb to blue and only compute the
        # rgb value for pixels that have been visited at least once
        # going to this trouble gave a ~2x speedup when hsv_to_rgb was 50% of the
        # overall compute time
        rgb = np.array([np.zeros(hue.shape), 
                        np.zeros(hue.shape), 
                        np.ones(hue.shape)*255]).transpose([2,1,0]).astype(int)

        # hsv[:,:,0] < 2/3 are the pixels with hue higher than 240/360=2/3, which
        # are the pixels that have been visited more than zero times
        rgb[hsv[:,:,0] < 2/3] = (hsv_to_rgb(hsv[hsv[:,:,0] < 2/3]) * 255).astype(int)

        # the values of imdata are ints packed with the r, g, and b channels
        imdata = (rgb[:,:,0] << 16) + (rgb[:,:,1] << 8) + rgb[:,:,2]
        """
        hue = (hue * 255).astype(int).transpose()
        imdata = self.packedColorLookup[hue]

        # draw a line down the middle to represent the meridian
        imdata[self.xCenter,:] = 0

        # outline the whole sky
        imdata[self.wholeSkyContour[:,0], self.wholeSkyContour[:,1]] = 0

        # draw a circle at zenith
        imdata[self.zenithAvoidContour[:,0], self.zenithAvoidContour[:,1]] = 0

        # draw the contour at 2 airmasses
        imdata[self.airmassContour[:,0], self.airmassContour[:,1]] = 0

        pygame.surfarray.blit_array(self.screen, imdata)
        pygame.display.flip()
        
        """ using matplotlib (which got ~4 frames/sec; pygame gets ~20 w/ imScale=2) 
        (also set contours to np.NaN not zero)
        # update the plotted image and colorbar
        self.image.set_data(imdata)
        self.image.set_clim(vmax=np.max(self.pixValues))

        plt.pause(0.001)
        #fig.canvas.draw()
        """

    def saveFrame(self, filename):
	pygame.image.save(self.screen, filename)



    #########################################################
    ## The methods below here are not used since I learned ##
    ## about _getChipNames or whatever it's called.        ##
    ## They could serve as an alternate potentially faster ##
    ## method though if the bugs were worked out.          ##
    #########################################################
    def getIncludedPixels(visit):
        ra = visit.ra
        dec = visit.dec
        rot = visit.rotation
        

        """
             ___________
            |           |
         ___|           |___
        |                   |
        |                   |
        |    o       i      |
        ||-------|C|----|   |
        |                   |
        |                   |
         ---             ---
            |           |
            |___________|

        """
        
        o = Telescope.fovWidth / 2
        i = Telescope.fovWidth * (3 / 5) / 2

        # these corners are [ra, dec] in radians
        bigSquare = np.array([
                       [ o,  o],
                       [ o, -o],
                       [-o, -o],
                       [-o,  o]
                    ])

        littleSquares = [
                np.array([[ o,  o], [ o,  i], [ i,  i], [ i,  o]]),
                np.array([[ o, -i], [ o, -o], [ i, -o], [ i, -i]]),
                np.array([[-i, -i], [-i, -o], [-o, -o], [-o, -i]]),
                np.array([[-i,  o], [-i,  i], [-o,  i], [-o,  o]]),
        ]

        allSquares = [bigSquare] + littleSquares

        # convert corners to [x, y, z]
        allSquares = [np.apply_along_axis(spherical2Cartesian, 1, square)
                      for square in allSquares]

        # rotate the FOV so it is centered on the visit (ra, dec)
        raAxis = np.array([0, 0, 1])
        allSquares = [np.array([rotate(corner, raAxis, ra) for corner in square])
                      for square in allSquares]

        decAxis = cross(spherical2Cartesian([ra, dec]), raAxis)
        allSquares = [np.array([rotate(corner, decAxis, dec) for corner in square])
                      for square in allSquares]

        # rotate the FOV about its center by visit.rotation
        rotAxis = spherical2Cartesian([ra, dec])
        allSquares = [np.array([rotate(corner, rotAxis, rot) for corner in square])
                      for square in allSquares]

        # convert back to (phi, theta)
        allSquares = [np.apply_along_axis(cartesian2Spherical, 1, square)
                      for square in allSquares]
        # project the corners using our existing wcs object
        # this returns [x,y] so we have to slice it with [:,[1,0]] to get [y,x]
        allSquares = [w.wcs_world2pix(square * RAD2DEG, 1)[:,[1,0]] * imScale
                      for square in allSquares]


        # generate candidate pixels which might fall in each square
        candidatePix = [getCandidatePix(square) for square in allSquares]

        # figure out which pixels are in each square
        includedPix = [[pix for pix in candidatePix[i] 
                        if isPixInSquare(pix, allSquares[i])]
                       for i in range(len(allSquares))]

        # store the included pixels in a set as tuples so we can do intersections
        includedPix = [set([tuple(pix) for pix in includedPix[i]])
                       for i in range(len(allSquares))]

        #print isPixInSquare([0,0],np.array([[1,1],[-1,1],[-1,-1],[1,-1]]))

        # return every pixel inside bigSquare that is in none of the little squares
        pixInLittleSquares = includedPix[1] | includedPix[2] | \
                             includedPix[3] | includedPix[4]

        #return pixInLittleSquares
        return includedPix[0] - pixInLittleSquares
        

    def getCandidatePix(corners):
        ys = corners[:,0]
        xs = corners[:,1]

        candidateXs = np.arange(round(min(xs)), round(max(xs)))
        candidateYs = np.arange(round(min(ys)), round(max(ys)))

        candidatePix = [[y, x] for x in candidateXs for y in candidateYs]
        return np.array(candidatePix)

    def rotate(v, axis, theta):
        """
        Code taken from user "unutbu"'s answer on stackoverflow: 
        """
        def rotation_matrix(axis, theta):
            """
            Return the rotation matrix associated with counterclockwise rotation about
            the given axis by theta radians.
            """
            axis = np.asarray(axis)
            theta = np.asarray(theta)
            axis = axis/np.sqrt(np.dot(axis, axis))
            a = np.cos(theta/2.0)
            b, c, d = -axis*np.sin(theta/2.0)
            aa, bb, cc, dd = a*a, b*b, c*c, d*d
            bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
            return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                             [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                             [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

        #def M(axis, theta):
        #	return expm3(cross(eye(3), axis/norm(axis)*theta))

        return dot(rotation_matrix(axis, theta), v)

    def isPixInSquare(point, square):
        """
        Crossing algorithm: run a semi-infinite array along the +x
        axis ending at pix and see how many times it crosses an edge
        of square. 

        Algorithm adapted from 
        http://www.realtimerendering.com/resources/GraphicsGems//gemsiv/ptpoly_haines/ptinpoly.c

        """

        # we're done once two edges have crossed the point's y coordinate
        # since a square (or any convex polygon) can only have that happen twice
        anyYCrossingsYet = False

        # every time the ray crosses an edge we negate this
        isInside = False

        # make indices more legible
        Y = 0
        X = 1

        for i in range(-1, len(square) - 1):
            # consider the square's edge joining vtx0 and vtx1
            vtx0 = square[i]
            vtx1 = square[i + 1]

            isY0AbovePoint = vtx0[Y] >= point[Y]
            isY1AbovePoint = vtx1[Y] >= point[Y]

            # only continue if the edge does not lie entirely
            # above or below point
            isPointBetweenYs = isY0AbovePoint != isY1AbovePoint
            if isPointBetweenYs:
                isX0RightOfPoint = vtx0[X] >= point[X]
                isX1RightOfPoint = vtx1[X] >= point[X]
                if isX0RightOfPoint == isX1RightOfPoint:
                    # if the edge lies entirely to the left or right of point
                    # then it intersects the ray exactly when it lies entirely
                    # to the right
                    if isX0RightOfPoint:
                        isInside = not isInside
                else:
                    # if the edge does not lie entirely to the left or right
                    # of point, we have to figure out whether the edge's
                    # intersection with the x axis lies to the left or right
                    # of point
                    xAxisIntersection = vtx1[X] - (vtx1[Y] - point[Y]) * \
                                        (vtx0[X] - vtx1[X]) / (vtx0[Y] - vtx1[Y])
                    if xAxisIntersection >= point[X]:
                        isInside = not isInside

                if anyYCrossingsYet:
                    return isInside
                anyYCrossingsYet = True

        return isInside


