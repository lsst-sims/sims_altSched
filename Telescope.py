from __future__ import division
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u

fovWidth = np.radians(3.5)
domSlitDiam = 2 * fovWidth
raftWidth = fovWidth / 5
minRotation = -np.pi/2
maxRotation = np.pi/2

# values from http://ops2.lsst.org/docs/current/system.html
maxAlt = np.radians(86.5)

# Kinematic and delay parameters for slew time computation

# speed in rads/sec
# acceleration in rads/second**2
domAltMaxSpeed = np.radians(1.75)
domAltAccel = np.radians(0.875)
domAltDecel = np.radians(0.875)

domAzMaxSpeed = np.radians(1.5)
domAzAccel = np.radians(0.75)
domAzDecel = np.radians(0.75)

telAltMaxSpeed = np.radians(3.5)
telAltAccel = np.radians(3.5)
telAltDecel = np.radians(3.5)
# assume accel == decel for calculations below
# (easy to change but they are the same anyway and I'm  lazy)
assert(telAltAccel == telAltDecel)

telAzMaxSpeed = np.radians(7.0)
telAzAccel = np.radians(7.0)
telAzDecel = np.radians(7.0)
assert(telAzAccel == telAzDecel)

# not used in slew calculation
Rotator_MaxSpeed = np.radians(3.5)
Rotator_Accel = np.radians(1.0)
Rotator_Decel = np.radians(1.0)

latitude = np.radians(-(30 + 14 / 60 + 40.7 / 3600))
longitude = np.radians(-(70 + 44 / 60 + 57.9 / 3600)) 

location = EarthLocation(lon = longitude * u.rad, lat = latitude * u.rad)

settleTime = 3

def calcSlewTime(altaz1, altaz2):
    # assume that we don't have to worry about the dome slew time
    # (might be reasonable if we never do long slews)
    # FYI this takes on the order of 10us for 1 slew calculation (longer now)

    # TODO also assumes we never max out the cable wrap-around constraint 
    deltaAlt = np.abs(altaz2[0] - altaz1[0])
    deltaAz  = np.abs(altaz2[1] - altaz1[1])

    deltaAz = min(deltaAz, np.abs(deltaAz - 2*np.pi))

    def uamSlewTime(d, vmax, a):
        # if you accelerate uniformely to telAltMaxSpeed
        # and then slow down uniformely to zero, you'll travel
        # a distance v_max^2 / a
        if d < vmax**2 / a:
            # to travel a distance d/2 while accelerating at a rate a,
            # it takes time sqrt(2(d/2)/a)
            slewTime = 2 * np.sqrt(d / a)
        else:
            # the time to accelerate/decelerate to/from v_max is 2v_max / a
            # and the distance covered in those two steps is v_max^2 / a
            # so the total time is the accel/decel time plus the remaining
            # distance over v_max
            slewTime = 2 * vmax / a + (d - vmax**2 / a) / vmax
        return slewTime

    #print "Delta alt", np.degrees(deltaAlt)
    #print "Delta az", np.degrees(deltaAz)
    telAltSlewTime = uamSlewTime(deltaAlt, telAltMaxSpeed, telAltAccel)
    telAzSlewTime  = uamSlewTime(deltaAz,  telAzMaxSpeed,  telAzAccel)

    # if we can fit both exposures in the dome slit, do so
    if deltaAlt**2 + deltaAz**2 < fovWidth**2:
        totDomTime = 0
    else:
        # else, we take the minimum time from two options:
        # 1. assume we line up alt in the center of the dome slit so we 
        #    minimize distance we have to travel in azimuth.
        # 2. line up az in the center of the slit
        # also assume:
        # * that we start out going maxspeed for both alt and az
        # * that we only just barely have to get the new field in the dome slit
        #   in one direction, but that we have to center the field in the other
        #   (which depends which of the two options used)
        # * that we don't have to slow down until after the shutter starts opening
        domDeltaAlt = deltaAlt

        # on each side, we can start out with the dome shifted away from the
        # center of the field by an amount domSlitRadius - fovRadius
        domDeltaAz = deltaAz - 2 * (domSlitDiam/2 - fovWidth/2)

        domAltSlewTime = domDeltaAlt / domAltMaxSpeed
        domAzSlewTime  = domDeltaAz  / domAzMaxSpeed

        totDomTime1 = max(domAltSlewTime, domAzSlewTime)

        domDeltaAlt = deltaAlt - 2 * (domSlitDiam/2 - fovWidth/2)
        domDeltaAz  = deltaAz
        domAltSlewTime = domDeltaAlt / domAltMaxSpeed
        domAzSlewTime  = domDeltaAz  / domAzMaxSpeed
        totDomTime2 = max(domAltSlewTime, domAzSlewTime)

        totDomTime = min(totDomTime1, totDomTime2)

    #print "slew alt/az", altSlewTime, azSlewTime
    #print
    return max(telAltSlewTime, telAzSlewTime, totDomTime) + settleTime
