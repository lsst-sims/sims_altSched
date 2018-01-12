from  __future__ import division
import numpy as  np
from lsst.sims.speedObservatory import Telescope
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from minis.RotationGenerator import RotationGenerator
from Visit import VisitPair
import config
import utils
from minis import Tiling

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
        tiling = Tiling.ThomsonTiling()
        allPointings = tiling.getTiling(config.tilingDensity)

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
