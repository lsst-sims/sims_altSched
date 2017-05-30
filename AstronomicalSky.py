from __future__ import division
import numpy as np
from astropy import units as u
import ephem
from datetime import datetime
from datetime import timedelta

import Config
from Telescope import Telescope

from matplotlib import pyplot as plt
import time


tel = ephem.Observer()
tel.lat = Telescope.latitude
tel.lon = Telescope.longitude
tel.horizon = "-12:00:00" # start observing when sun is 12 deg below horizon
sun = ephem.Sun()

starts = {}
ends = {}

def nightLength(nightNum):
    return nightEnd(nightNum) - nightStart(nightNum)

def nightStart(nightNum):
    if nightNum in starts:
        return starts[nightNum]
    startDate = datetime.fromtimestamp(Config.surveyStartTime)
    startDate = startDate.replace(hour=12, minute=0, second=0)
    curDate = startDate + timedelta(nightNum)
    dublinJD = tel.next_setting(sun, start=curDate)
    mjd = dublinJD + 15019.5 # convert from dublin JD to modified JD
    unix = (mjd - 40587) * 86400
    starts[nightNum] = unix
    return unix

def nightEnd(nightNum):
    if nightNum in ends:
        return ends[nightNum]
    startDate = datetime.fromtimestamp(Config.surveyStartTime)
    startDate = startDate.replace(hour=12, minute=0, second=0)
    curDate = startDate + timedelta(nightNum)
    dublinJD = tel.next_rising(sun, start=curDate)
    mjd = dublinJD + 15019.5 # convert from dublin JD to modified JD
    unix = (mjd - 40587) * 86400
    ends[nightNum] = unix
    return unix

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


def localSiderialTime(longitude, times):
    # this is the timestamp of 12pm 1/1/2000 
    J2000 = 946728000
    daysSinceJ2000 = (times - J2000) / 3600 / 24
    dts = np.vectorize(datetime.utcfromtimestamp)(times)
    def getTimeInHours(dt):
        return dt.hour + dt.minute / 60 + dt.second / 3600
    timesInHours = np.vectorize(getTimeInHours)(dts)

    localSiderialTimes = 100.46 + 0.985647 * daysSinceJ2000 + \
                        np.degrees(longitude) + 15 * timesInHours
    localSiderialTimes %= 360
    localSiderialTimes = np.radians(localSiderialTimes)
    return localSiderialTimes

# precompute HA/dec => alt/az lookup table
# TODO I expected this lookup table to make radec2altaz much faster, but it
# only made it marginally faster, if at all. Should look into this with
# a line profiler at some point
pi = np.pi
sin = np.sin
cos = np.cos

# at this resolution it takes ~0.8sec to compute this lookup table
# and for points far from singular points in the map, the deviation
# between adjacent lookup table pixels is <~0.1 degrees
g_nHAs = 3200
g_nDecs = 1600

g_HA, g_dec = np.mgrid[0     : 2*pi : 2*pi/g_nHAs,
                       -pi/2 : pi/2 : pi/g_nDecs]

g_sinDec = sin(g_dec)
g_sinAlt = g_sinDec * sin(Telescope.latitude) + \
         cos(g_dec) * cos(Telescope.latitude) * cos(g_HA)
g_alt = np.arcsin(g_sinAlt)

g_cosA = (g_sinDec - g_sinAlt * sin(Telescope.latitude)) / \
               (cos(g_alt) * cos(Telescope.latitude))
g_cosA = np.clip(g_cosA, -1, 1)
g_a = np.arccos(g_cosA)

g_sinHA = sin(g_HA)
g_az = g_a
g_az[g_sinHA >= 0] = 2*pi - g_a[g_sinHA >= 0]

g_altLookup = g_alt
g_azLookup = g_az

"""
plt.title("alt")
plt.imshow(g_altLookup)
plt.figure()
plt.title("az")
plt.imshow(g_azLookup)
plt.show()
"""

def radec2altaz(radec, times):
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
    astropyResults = np.vstack([a.alt.rad, a.az.rad]).T
    #return astropyResults
    """

    # the method below returns values within ~1 degree of what the astropy
    # code gives (for reasonable dates), and it takes ~30us instead of ~15ms
    # TODO sometimes they differ by >2 degrees in the north -- this is too much
    # maybe interpolate? or use higher res lookup?

    # formulas from http://www.stargazing.net/kepler/altaz.html
    
    LSTs = localSiderialTime(Telescope.longitude, times)
    if not isinstance(LSTs, np.ndarray):
        LSTs = LSTs * np.ones(ra.shape)
    HA = (LSTs - ra) % (2*np.pi)

    HAIndex = (HA / (2*np.pi) * g_nHAs).astype(int)
    decIndex = ((dec+np.pi/2) / np.pi * g_nDecs).astype(int)

    alt = g_altLookup[HAIndex, decIndex]
    az = g_azLookup[HAIndex, decIndex]

    #diff = np.degrees(astropyResults - np.vstack([alt, az]).T)
    #print "diff", diff[diff > 0.5]
    return np.vstack([alt, az]).T


def altaz2radec(altaz, times):
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
    
    LSTs = localSiderialTime(Telescope.longitude, times)
    sin = np.sin
    cos = np.cos
    lat = Telescope.latitude
    sinDec = sin(alt) * sin(lat) + cos(alt) * cos(lat) * cos(az)
    dec = np.arcsin(sinDec)
    sinHourAngle = -1 * sin(az) * cos(alt) / cos(dec)
    cosHourAngle = (sin(alt) - sin(dec) * sin(lat)) / (cos(dec) * cos(lat))
    hourAngle = np.arcsin(sinHourAngle)
    hourAngle[cosHourAngle <= 0] = np.pi - hourAngle[cosHourAngle < 0]
    ra = LSTs - hourAngle

    # this method runs in ~40us vs ~12ms for the astropy code and returns
    # results within ~1deg as above
    #print "diff", np.degrees(astropyResults) - np.degrees(np.vstack([ra, dec]).T)
    return np.vstack([ra, dec]).T
