from __future__ import division

import time

import numpy as np
import ephem
from datetime import datetime
from datetime import timedelta

import Config
from Telescope import Telescope
import Utils

from matplotlib import pyplot as plt
import palpy


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
    #return Config.surveyStartTime + nightNum * 3600*24
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
    #return Config.surveyStartTime + nightNum * 3600*24 + 12 * 3600
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
    altaz = np.array([[np.pi/2, 0]])
    radec = altaz2radec(altaz, time)
    return radec[0,0]

def radecOfMoon(time):
    moon = ephem.Moon(Utils.mjd2djd(Utils.unix2mjd(time)))
    return (moon.ra, moon.dec)

def phaseOfMoon(time):
    moon = ephem.Moon(Utils.mjd2djd(Utils.unix2mjd(time)))
    return moon.moon_phase

def getExpTime(ra, dec, otherstuff = None):
    return 30


def unix2lst(longitude, time):
    mjd = Utils.unix2mjd(time)
    lst = palpy.gmst(mjd) + longitude
    lst %= 2*np.pi
    return lst

def radec2altaz(radec, time):
    # code adapted from lsst-ts/ts_astrosky_model and lsst-ts/ts_dateloc
    ra = radec[:,0]
    dec = radec[:,1]
    lst = unix2lst(Telescope.longitude, time)
    ha = lst * np.ones(ra.shape) - ra
    az, alt = palpy.de2hVector(ha, dec, Telescope.latitude)
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
    # TODO do this with palpy
    
    lst = unix2lst(Telescope.longitude, time)

    sin = np.sin
    cos = np.cos
    lat = Telescope.latitude
    sinDec = sin(alt) * sin(lat) + cos(alt) * cos(lat) * cos(az)
    dec = np.arcsin(sinDec)
    sinHourAngle = -1 * sin(az) * cos(alt) / cos(dec)
    cosHourAngle = (sin(alt) - sin(dec) * sin(lat)) / (cos(dec) * cos(lat))
    hourAngle = np.arcsin(sinHourAngle)
    hourAngle[cosHourAngle <= 0] = np.pi - hourAngle[cosHourAngle < 0]
    ra = lst - hourAngle

    #print "diff", np.degrees(astropyResults) - np.degrees(np.vstack([ra, dec]).T)
    return np.vstack([ra, dec]).T
