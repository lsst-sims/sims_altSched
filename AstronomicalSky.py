from __future__ import division
import numpy as np
import Telescope
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation
from astropy.coordinates import ICRS
from astropy import units as u
from datetime import datetime

import Config

def nightLength(nightNum):
    return 60 * 60 * 8

def nightStart(nightNum):
    return Config.surveyStartTime + nightNum * 60 * 60 * 24

def nightEnd(nightNum):
    return nightStart(nightNum) + nightLength(nightNum)

def nightNum(time):
    return int((time - Config.surveyStartTime) / 3600 / 24)

def raOfMeridian(time):
    """
    oneDay = 60 * 60 * 24
    oneYear = oneDay * 365.25

    dayContribution = (time % oneDay) * 2*np.pi / oneDay 
    yearContribution = (time % oneYear) * 2*np.pi / oneYear
    return (dayContribution + yearContribution) % (2*np.pi)
    """
    t = datetime.utcfromtimestamp(time)
    altaz = np.array([[np.pi/2, 0]])
    radec = altaz2radec(altaz, time)
    return radec[0,0]

def getExpTime(ra, dec, otherstuff = None):
    return 30


def localSiderialTime(longitude, time):
    # this is the timestamp of 12pm 1/1/2000 
    J2000 = 946728000
    daysSinceJ2000 = (time - J2000) / 3600 / 24
    dt = datetime.utcfromtimestamp(time)
    timeInHours = dt.hour + dt.minute / 60 + dt.second / 3600

    localSiderialTime = 100.46 + 0.985647 * daysSinceJ2000 + \
                        np.degrees(longitude) + 15 * timeInHours
    localSiderialTime %= 360 
    localSiderialTime = np.radians(localSiderialTime)
    return localSiderialTime

def radec2altaz(radec, time):
    ra = radec[:,0]
    dec = radec[:,1]

    """
    This is unbelievably inexplicably slow (>~90% runtime), ~15ms/call
    I made a call graph for the program using this method and it called
    a zillion random functions none of which took a significant amount of 
    time individually...

    t = datetime.utcfromtimestamp(time)
    altazframe = AltAz(location = Telescope.location, obstime = t)

    icrs = ICRS(ra=ra * u.rad, dec=dec * u.rad)
    a = icrs.transform_to(altazframe)
    astropyResults = np.degrees(np.vstack([a.alt.rad, a.az.rad]).T)
    """

    # the method below returns values within ~1 degree of what the astropy
    # code gives (for reasonable dates), and it takes ~30us instead of ~15ms

    # formulas from http://www.stargazing.net/kepler/altaz.html
    
    LST = localSiderialTime(Telescope.longitude, time)
    hourAngle = (LST * np.ones(ra.shape) - ra) % (2*np.pi)
    sin = np.sin
    cos = np.cos
    sinalt = sin(dec) * sin(Telescope.latitude) + \
             cos(dec) * cos(Telescope.latitude) * cos(hourAngle)
    alt = np.arcsin(sinalt)

    cosA = (sin(dec) - sin(alt) * sin(Telescope.latitude)) / \
                      (cos(alt) * cos(Telescope.latitude))
    a = np.arccos(cosA)
    az = a
    az[sin(hourAngle) >= 0] = 2*np.pi - a[sin(hourAngle) >= 0]

    #print "diff", astropyResults - np.degrees(np.vstack([alt, az]).T)
    return np.vstack([alt, az]).T


def altaz2radec(altaz, time):
    alt = altaz[:,0]
    az = altaz[:,1]

    """
    #Same problem as radec2altaz...
    t = datetime.utcfromtimestamp(time)
    altaz = AltAz(location = Telescope.location, obstime = t, alt = alt * u.rad, az = az * u.rad)
    radec = altaz.transform_to(ICRS)
    astropyResults = np.vstack([radec.ra.rad, radec.dec.rad]).T
    """

    # formulas from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
    
    LST = localSiderialTime(Telescope.longitude, time)
    sin = np.sin
    cos = np.cos
    lat = Telescope.latitude
    sinDec = sin(alt) * sin(lat) + cos(alt) * cos(lat) * cos(az)
    dec = np.arcsin(sinDec)
    sinHourAngle = -1 * sin(az) * cos(alt) / cos(dec)
    cosHourAngle = (sin(alt) - sin(dec) * sin(lat)) / (cos(dec) * cos(lat))
    hourAngle = np.arcsin(sinHourAngle)
    hourAngle[cosHourAngle <= 0] = np.pi - hourAngle[cosHourAngle < 0]
    ra = LST - hourAngle

    # this method runs in ~40us vs ~12ms for the astropy code and returns
    # results within ~1deg as above
    #print "diff", np.degrees(astropyResults) - np.degrees(np.vstack([ra, dec]).T)
    return np.vstack([ra, dec]).T
