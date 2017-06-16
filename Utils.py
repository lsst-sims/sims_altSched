from __future__ import division

import numpy as np
from astropy.time import Time
import multiprocessing as mp
from multiprocessing import Pool
import itertools

def areRasInRange(ras, (minRa, maxRa)):
    if minRa < maxRa:
        return (minRa < ras) & (ras < maxRa)
    else:
        return (minRa < ras) | (ras < maxRa)

def isRaInRange(ra, (minRa, maxRa)):
    if minRa < maxRa:
        return minRa < ra and ra < maxRa
    else:
        return minRa < ra or  ra < maxRa

def unix2mjd(timestamp):
    return timestamp / 86400 + 40587

def mjd2unix(mjd):
    return (mjd - 40587) * 86400

def mjd2djd(mjd):
    """
    Convert Modified Julian Date to Dublin Julian Date (what pyephem uses).
    (this function adapted from Peter Yoachim's code)
    """
    doff = 15019.5 # this equals ephem.Date(0)-ephem.Date('1858/11/17')
    return mjd - doff

def timetsamp(mjd):
    t = Time(mjd, format="mjd")
    return t.unix

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


# thought I would need this but seems like
# list(itertools.chain.from_iterable(array)) is fast enough
# by itself. Leaving this here just in case...
def _sumWorker(array):
    return list(itertools.chain.from_iterable(array))

def sumArrayOfLists(arrayOfLists):
    # 1000 arbitrarily chosen
    sublistSize = 1000

    nCores = mp.cpu_count()
    p = Pool(nCores)
    nSublists = int(len(arrayOfLists) / sublistSize)
    summed = p.map(_sumWorker, np.array_split(arrayOfLists, nSublists))
    p.close()
    p.join()

    return _sumWorker(summed)
