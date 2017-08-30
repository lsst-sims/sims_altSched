from __future__ import division

import numpy as np
from astropy.time import Time
import multiprocessing as mp
from multiprocessing import Pool
import itertools
from lsst.sims.speedObservatory import Telescope
import Config
from Config import NORTH, SOUTH, EAST, SOUTHEAST

def areRasInRange(ras, raRange):
    """ Calculate which of `ras` are within `raRange`

    Parameters
    ----------
    ras : numpy.ndarray
        Array of RAs to be tested
    raRange : 2-tuple of (minRa, maxRa)
        The range to be tested

    Returns
    -------
    A numpy.ndarray indicating which of `ras` are within `raRange`
    """
    minRa, maxRa = raRange
    if minRa < maxRa:
        return (minRa < ras) & (ras < maxRa)
    else:
        return (minRa < ras) | (ras < maxRa)

def isRaInRange(ra, raRange):
    """ Calculate whether `ra` is within `raRange`

    Parameters
    ----------
    ra : float
        The RA to be tested
    raRange : 2-tuple of (minRa, maxRa)
        The range to be tested

    Returns
    -------
    A boolean indicating whether `ra` is within `raRange`
    """

    minRa, maxRa = raRange
    if minRa < maxRa:
        return minRa < ra and ra < maxRa
    else:
        return minRa < ra or  ra < maxRa

def directionOfDec(dec):
    """ Calculate which cardinal direction `dec` is in

    Parameters
    ----------
    dec : float
        The declination of interest in radians

    Returns
    -------
    One of NORTH, SOUTH, or EAST, where EAST means the dec
    falls within the zenith dec band.

    Notes
    -----
    Calling things in the zenith dec band EAST is kind of stupid, sorry
    """
    if dec > Config.maxDec or dec < Config.minDec:
        raise ValueError("Provided dec of " + str(dec) + " is outside " + \
                         "of the survey area.")
    if dec > Telescope.latitude + Config.zenithBuffer:
        return NORTH
    elif dec > Telescope.latitude - Config.zenithBuffer:
        return EAST
    else:
        return SOUTH

def areaInDir(direction):
    """ Returns how much area (in steradians) falls in `direction

    Parameters
    ----------
    direction : enum(NORTH, SOUTH, EAST, SOUTHEAST)
        The direction of interest

    Returns
    -------
    The amount of survey area in `direction`, where EAST means
    "within the zenith dec band" (see note from directionOfDec())
    """
    buf = Config.zenithBuffer
    # int_{dec1}^{dec2} 2\pi r dr, where r=\cos\theta
    # = 2\pi (\sin(dec2) - \sin(dec1))
    if direction == NORTH:
        return 2*np.pi*(np.sin(Config.maxDec) -
                        np.sin(Telescope.latitude + buf))

    elif direction == SOUTH:
        return 2*np.pi*(np.sin(Telescope.latitude - buf) -
                        np.sin(Config.minDec))
    elif direction == EAST:
        return 2*np.pi*(np.sin(Telescope.latitude + buf) -
                        np.sin(Telescope.latitude - buf))

    elif direction == SOUTHEAST:
        return areaInDir(SOUTH) + areaInDir(EAST)
    else:
        raise ValueError("Invalid direction " + str(direction))

# I could use one of astropy's implementation of this
# but one of them is 10x slower than this and the 
# other requires using Longitude or Angle classes or
# something and I don't want to go down that road...
def spherical2Cartesian(phi, theta):
    x = np.cos(theta) * np.cos(phi)
    y = np.cos(theta) * np.sin(phi)
    z = np.sin(theta)

    if isinstance(z, np.ndarray):
        return np.vstack([x, y, z]).T
    else:
        return [x, y, z]
