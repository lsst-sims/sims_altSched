from  __future__ import division
import numpy as  np
from lsst.sims.speedObservatory import Telescope
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from minis.RotationGenerator import RotationGenerator
from Visit import VisitPair
import config
import utils
from astropy import wcs

class MiniSurvey:
    """ This is a static class that generates tilings

    The main public method is newMiniSurvey, which returns a set of
    VisitPairs in a requested RA/Dec rectangle.

    TODO a refactoring I've been meaning to do is to make a base Tiling
    interface with a method `getTiling(density)` that returns a full-sky
    tiling with the specified number of pointings/steradian. Then subclasses
    could execute different tilings.

    It's important to be able to choose the density, since if the density of
    a tiling is too high, the scheduler won't be able to keep up with the
    meridian or will have to skip visits, which would likely increase slew
    times. If the density is too low, the scheduler will either get ahead
    of the meridian or will have to double cover, which wastes time observing
    a field more than twice per night.
    """

    rotationGenerator = RotationGenerator()

    def __init__(self):
        raise NotImplementedError("MiniSurvey is static")

    @classmethod
    def setLatestRotation(cls, rotation, direction):
        # called when resuming the survey from a checkpoint
        # TODO but checkpointing is not implemented
        cls.rotationGenerator = RotationGenerator(rotation, direction)

    @classmethod
    def newMiniSurvey(cls, telescope, minRa, maxRa, direction):
        """ Return a new tiling of the sky

        Parameters
        ----------
        telescope : lsst.sims.speedObservatory.Telescope
            The telescope that the minisurvey is for
        minRa : float
            The minimum RA that should be included in the minisurvey (radians)
        maxRa : float
            The maximum RA that should be included in the minisurvey (radians)
        direction : enum(NORTH, SOUTH, EAST)
            The cardinal direction the minisurvey should be in

        Returns
        -------
        A set of visitpairs that tile the requested area

        Notes
        -----
        min/maxRa should be the RA range for the night. This method will
        shift over the RA range to the right if direction is EAST. I know this
        is non-ideal, sorry.
        """

        # TODO rotation isn't fully implemented yet
        rotation = next(cls.rotationGenerator.rotations(telescope))

        # get a tiling over the whole sphere
        allPointings = cls._realGeneratePointings(telescope)

        # now take the pointings and rotate them around the z and x
        # axes so they are randomly dithered on the sphere
        zDither = np.random.randn()
        xDither = np.random.randn()
        
        # rotating around z is easy
        allPointings[:,0] += zDither

        # rotating around x is harder
        sin = np.sin
        cos = np.cos
        phi = allPointings[:,0]
        theta = allPointings[:,1]
        newPhi = np.arctan2(-sin(xDither)*sin(theta) + cos(xDither)*cos(theta)*sin(phi),
                            cos(theta)*cos(phi))
        newPhi %= 2 * np.pi
        newTheta = np.arcsin(cos(xDither)*sin(theta) + sin(xDither)*cos(theta)*sin(phi))
        allPointings = np.vstack([newPhi, newTheta]).T

        # TODO should rotate around y as well or else you get structured
        # noise in the final result

        # calculate min/maxDec that correspond to the passed-in direction
        # and modify min/maxRa if we're in the zenith dec band
        if direction == config.NORTH:
            minDec = telescope.latitude + config.zenithBuffer
            maxDec = config.maxDec
            # min/maxRa remain unchanged
        elif direction == config.SOUTH:
            minDec = config.minDec
            maxDec = telescope.latitude - config.zenithBuffer
            # min/maxRa remain unchanged
        elif direction == config.EAST:
            minDec = telescope.latitude - config.zenithBuffer
            maxDec = telescope.latitude + config.zenithBuffer
            minRa += config.zenithBuffer + config.zenithBufferOffset
            maxRa += config.zenithBuffer + config.zenithBufferOffset
        else:
            raise ValueError("Invalid direction: " + str(direction))

        # choose the subset of pointings that lie in the min/maxRa/Dec rectangle
        validRa = utils.areRasInRange(allPointings[:,0], (minRa, maxRa))
        validDec = ((minDec < allPointings[:,1]) &
                    (maxDec > allPointings[:,1]))
        validMask = validRa & validDec

        pointings = allPointings[np.where(validMask)]

        # create visitpairs from the calculated pointings
        visitPairs = [VisitPair(pointing[0], pointing[1], rotation)
                      for pointing in pointings]

        return set(visitPairs)

    @classmethod
    def _generateRandomPointings(self, telescope):
        numPointings = int(4 * np.pi / (telescope.fovWidth**2))
        thetas = np.arcsin(np.random.random(numPointings)*2 - 1)
        phis = np.random.random(numPointings)*2*np.pi

        pointings = np.vstack([thetas, phis]).transpose()
        return pointings


    @classmethod
    def _realGeneratePointings(self, telescope):
        """ Tiles the sky with a 2-raft-overlap pattern

        Or rather, the tiling has 2 rafts of overlap in flat space. Turns
        out that after projecting onto the sphere, this isn't really the
        case any more. Fortunately, since we do large-scale full-tiling
        dithering, the exact tiling doesn't matter all that much.
        """
        # to generate the pointings, start with the pointing pattern in flat
        # space, then map the pattern onto a quadrilateralized spherical cube

        # the coordinates in flat space are (u,v) where u and v represent
        # the number of rafts of displacement in the horizontal and vertical
        # directions, respectively, from the bottom left corner (0,0)

        # in these coordinates, focal plane is a square with side length
        # equal to 5 (since the focal plane is 5x5 rafts)

        # although the raft and fov widths are not hard-coded, the number of 
        # rafts is hard-coded. After all, this algorithm completely relies on
        # the fact that the focal plane is a 5x5 square with 1x1 squares
        # removed from the corners

        # I'm guessing there's a better way to do this, but here goes...

        ### TODO XXX IMPORTANT XXX TODO ###
        ### multiplying `raftWidth` by some float near 1 is a hacky way
        ### to change the spacing of the tiling. A way that allows
        ### the user to choose, say, the number of pointings/steradian
        ### would be better but this works
        ### TODO XXX IMPORTANT XXX TODO ###
        raftWidth = telescope.raftWidth * 0.98

        cubeSideLength = np.arccos(1/3) # in radians
        numRaftsPerCubeSide = int(cubeSideLength / raftWidth)

        minU = 0
        minV = 0

        # we need to tile enough to be sure that we can completely cover
        # the unfolded cube with the pattern
        #maxU = int(2 * np.pi / raftWidth)
        #maxV = int(np.pi / raftWidth)

        # if you unfold the cube like this:
        #  _
        # | |_____
        # |  _____|
        # |_|
        # then we need pattern a rectangle of size numRaftsPerCubeSide * (4,3)
        # I'm using (5,4) instead to avoid the edges of the pattern
        maxU = 5 * numRaftsPerCubeSide 
        maxV = 4 * numRaftsPerCubeSide

        u = minU
        v = minV
        pointings = []

        # this while loop populates the pointings array with (u,v) tuples 
        # in the proper pointing pattern. The pattern fills a rectangle
        # stretching from (minU, minV) to (maxU, maxV)
        while True:
            prevBottomV = v
            prevBottomU = u
            while True:
                # scan up and to the left, moving left 1 raft and up 4 rafts 
                # each time

                pointings.append((u,v))

                u -= 1
                v += 4 

                if u < minU or v > maxV:
                    # we're done with this column
                    break
                    
            if prevBottomV >= 1: 
                # start the next column shifted 5 rafts right and one raft down
                u = prevBottomU + 5 
                v = prevBottomV - 1
            else:
                # start the next column shifted 4 rafts right and 3 rafts up
                u = prevBottomU + 4
                v = prevBottomV + 3

            # now check if we've reached the bottom right corner
            # if so, we need to start the next column farther up
            if u > maxU:
                v += 4 * (u - maxU)
                u = maxU

                # if v is now larger than maxV, that means we're done
                if v > maxV:
                    break

        pointings = np.array(pointings).astype("float64")

        
        # astropy lays out the cube like this
        # _
        #|0|_____
        #|1 4_3_2|
        #|5|
        # 
        # where (0,0) is (annoyingly) in the center of face 1
        # and where faces 4, 3, and 2 are repeated to the left of
        # face 1
        #     

        # convert to "degrees." Each cube face has a side length of 90 degrees 
        pointings *= 90 / numRaftsPerCubeSide

        # shift the origin so the pattern covers the entire unfolded cube
        # (and doesn't double cover anything)
        pointings[:,0] -= 45
        pointings[:,1] -= 135 

        # convert to RA/DEC
        w = wcs.WCS(naxis=2)
        # QSC is quadrilateralized spherical cube
        w.wcs.ctype = ["RA---QSC", "DEC--QSC"]
        # the 2nd argument is the "origin" of the coordinate system
        # "1" represents (1,1) and I think it conveys that the QSC 
        # point (0,0) corresponds to RA/DEC of (0,0)? TODO
        pointings = w.wcs_pix2world(pointings, 0) 
        
        # nan indicates the pointing was not on any cube face
        pointings = pointings[~np.isnan(pointings[:,0])]

        # and convert back to radians
        pointings = np.radians(pointings)

        return pointings
