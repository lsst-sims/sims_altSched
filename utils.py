from __future__ import division

import numpy as np
from astropy.time import Time
import multiprocessing as mp
from multiprocessing import Pool
import itertools
from lsst.sims.speedObservatory import Telescope
import config
from config import NORTH, SOUTH, EAST, SOUTHEAST
from lsst.sims.speedObservatory import sky

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
    if dec > config.maxDec or dec < config.minDec:
        raise ValueError("Provided dec of " + str(dec) + " is outside " + \
                         "of the survey area.")
    if dec > Telescope.latitude + config.zenithBuffer:
        return NORTH
    elif dec > Telescope.latitude - config.zenithBuffer:
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
    buf = config.zenithBuffer
    # int_{dec1}^{dec2} 2\pi r dr, where r=\cos\theta
    # = 2\pi (\sin(dec2) - \sin(dec1))
    if direction == NORTH:
        return 2*np.pi*(np.sin(config.maxDec) -
                        np.sin(Telescope.latitude + buf))

    elif direction == SOUTH:
        return 2*np.pi*(np.sin(Telescope.latitude - buf) -
                        np.sin(config.minDec))
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

def cartesian2Spherical(x, y, z):
    phi = np.arctan2(y, x) + np.pi
    theta = np.arctan2(np.sqrt(x**2+y**2), z) - np.pi/2

    if isinstance(z, np.ndarray):
        return np.vstack([phi, theta]).T
    else:
        return [phi, theta]


def _getTimeHorizon(ra,dec,limit,tmin,tmax,limit_type='lower',n_inter=10):
    """ Estimate times for which the field is visible (above or below a given line)
          Input : ra,dec of the field
                       lim: the limit below or above which the field is
                       tmin,tmax= min and max times of observation  
                       limit_type : below ('upper') or above ('lower') the limit
          Output: times of min and max when the field crosses the limit
    """
    tmin_h = tmin
    tmax_h = tmax

    while (tmax_h - tmin_h) > 10.:

        times_h = np.linspace(tmin_h, tmax_h, n_inter)
        cross = False
        for i in range(len(times_h) - 1):
            alt_a, az_a = sky.radec2altaz(ra, dec, times_h[i])
            alt_b, az_b = sky.radec2altaz(ra, dec, times_h[i+1])

            diff_a=alt_a - limit
            diff_b=alt_b - limit

            if (diff_a*diff_b) <= 0:
                tmin_h = times_h[i]
                tmax_h = times_h[i+1]
                cross = True
                break

        if cross is False:
            #  No horizon line crossed -> get out !
            alt, az = sky.radec2altaz(ra, dec, times_h[0])
            if (((alt - limit) >= 0 and limit_type == 'lower') or
                ((alt - limit) <= 0 and limit_type == 'upper')):
                return (np.min(times_h), np.max(times_h))
            else:
                return (-1.0, -1.0)
            break

    alt_a, az_a = sky.radec2altaz(ra, dec,tmin_h)
    alt_b, az_b = sky.radec2altaz(ra, dec,tmax_h)
    if  limit_type == 'lower':
        if (alt_a - limit)>=0:
            return (tmin, tmin_h)
        else:
            return (tmax_h, tmax)
    else:
        if (alt_a-limit) >= 0:
            return (tmin_h, tmax)
        else:
            return (tmin, tmax_h)

def _getMeridianMin(ra, dec, tmin, tmax, telescope, n_inter=100):
    """ Estimate times for which the field is closest to the meridian
          Input : ra,dec of the field
                       tmin,tmax= min and max times of observation
                       telescope : the telescope
          Output: tmin, tmax: times of min and max of observation
                  ha: minimal hour angle  of the field (between tmin and tmax)
                  time_fi : time corresponding to the minimun ha
                            (between tmin and tmax)
    """
    tmin_h = tmin
    tmax_h = tmax

    while (tmax_h - tmin_h) > 10.:

        times_h = np.linspace(tmin_h, tmax_h, n_inter)
        mini = False
        for i in range(len(times_h) - 1):
            lst_a = sky.unix2lst(telescope.longitude, times_h[i])
            ha_a = lst_a - ra
            lst_b = sky.unix2lst(telescope.longitude, times_h[i + 1])
            ha_b = lst_b - ra
            if ha_a > np.pi:
                ha_a = 2.*np.pi-ha_a
            if ha_b > np.pi:
                ha_b = 2.*np.pi-ha_b
            
            if (ha_a - ha_b) <=0:
                tmin_h = times_h[i]
                tmax_h = times_h[i+1]
                mini = True
                break

        if not mini:
            lst = sky.unix2lst(telescope.longitude, times_h[len(times_h) - 1])
            ha = lst - ra
            if ha > np.pi:
                ha = 2.*np.pi-ha
            return tmin, tmax, ha, times_h[len(times_h) - 1]

    time_fi = np.mean([tmin_h, tmax_h])
    lst = sky.unix2lst(telescope.longitude, time_fi)
    ha = np.abs(lst - ra)
    if ha > np.pi:
        ha=2.*np.pi - ha
    return tmin, tmax, ha, time_fi

def getTimeObs(ra,dec,exptime,telescope,tmin,tmax):
    """ Estimate times for which the field is visible (above the horizon)
          Input: ra,dec : ra,dec of the field
                 exptime : total exposure time of the field
                 telescope : telescope model
                 tmin, tmax : min and max times of observation
         Output: tmin, tmax : min and max times for the field being above
                              horizon and completely observable (ie the
                              complete sequence of observation can be done)
                        ha : min hour angle between tmin and tmax
                        time_fi : time of min(ha)
    """
    tmin_up, tmax_up = _getTimeHorizon(ra, dec, telescope.minAlt, tmin, tmax)
    tmin_down, tmax_down = _getTimeHorizon(ra, dec, telescope.maxAlt,
                                           tmin, tmax, limit_type='upper')

    if tmin_up * tmax_up == 1 or  tmin_down * tmax_down == 1:
        return -1.0, -1.0, -1.0, -1.0

    if tmin_down > tmax_up or tmin_up > tmax_down:
        # No recovery between up and down region
        return -1.0, -1.0, -1.0, -1.0
    else:
        tmin_fi = tmin_down
        tmax_fi = min(tmax_up, tmax_down)
        if tmin_up >= tmin_down:
            tmin_fi = tmin_up
        #check whether we will have time to perform the full sequence
        if tmax_fi - tmin_fi < exptime:
            return tmin_fi, tmax_fi, -1.0, -1.0
        #we will not have the time to perform the sequence before the end of the night
        if tmax-tmin_fi < exptime:
            return tmin_fi, tmax_fi, -1.0, -1.0

        return _getMeridianMin(ra, dec, tmin_fi, tmax_fi-exptime, telescope)
