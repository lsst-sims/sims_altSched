import numpy as np
from astropy.time import Time

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

def mjd(timestamp):
    t = Time(timestamp, format="unix")
    return t.mjd

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
