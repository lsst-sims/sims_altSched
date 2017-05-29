from __future__ import division
import numpy as np
from astropy.coordinates import EarthLocation
from astropy import units as u

class Telescope:
    # non-configurable (AstronomicalSky uses these and I don't
    # want to slow it down by requiring it create a Telescope object)
    # (don't modify these on an instance)
    latitude = np.radians(-(30 + 14 / 60 + 40.7 / 3600))
    longitude = np.radians(-(70 + 44 / 60 + 57.9 / 3600))
    filters = ["u", "g", "r", "i", "z", "y"]
    filterId = {}
    for i, filter in enumerate(filters):
        filterId[filter] = i

    def __init__(self):
        self.fovWidth = np.radians(3.5)
        self.domSlitDiam = 2 * self.fovWidth
        self.raftWidth = self.fovWidth / 5
        self.minRotation = -np.pi/2
        self.maxRotation = np.pi/2

        # values from http://ops2.lsst.org/docs/current/system.html
        self.maxAlt = np.radians(86.5)

        # Kinematic and delay parameters for slew time computation

        # speed in rads/sec
        # acceleration in rads/second**2
        self.domAltMaxSpeed = np.radians(1.75)
        self.domAltAccel = np.radians(0.875)
        self.domAltDecel = np.radians(0.875)

        self.domAzMaxSpeed = np.radians(1.5)
        self.domAzAccel = np.radians(0.75)
        self.domAzDecel = np.radians(0.75)

        self.telAltMaxSpeed = np.radians(3.5)
        self.telAltAccel = np.radians(3.5)
        self.telAltDecel = np.radians(3.5)
        # assume accel == decel for calculations below
        # (easy to change but they are the same anyway and I'm  lazy)
        assert(self.telAltAccel == self.telAltDecel)

        self.telAzMaxSpeed = np.radians(7.0)
        self.telAzAccel = np.radians(7.0)
        self.telAzDecel = np.radians(7.0)
        assert(self.telAzAccel == self.telAzDecel)

        # not used in slew calculation
        self.Rotator_MaxSpeed = np.radians(3.5)
        self.Rotator_Accel = np.radians(1.0)
        self.Rotator_Decel = np.radians(1.0)

        self.location = EarthLocation(lon = self.longitude * u.rad,
                                      lat = self.latitude * u.rad)

        self.settleTime = 3
        self.filterChangeTime = 120

    def calcSlewTime(self, altaz1, filter1, altaz2, filter2):
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
        telAltSlewTime = uamSlewTime(deltaAlt, self.telAltMaxSpeed, self.telAltAccel)
        telAzSlewTime  = uamSlewTime(deltaAz,  self.telAzMaxSpeed,  self.telAzAccel)

        # if we can fit both exposures in the dome slit, do so
        if deltaAlt**2 + deltaAz**2 < self.fovWidth**2:
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
            domDeltaAz = deltaAz - 2 * (self.domSlitDiam/2 - self.fovWidth/2)

            domAltSlewTime = domDeltaAlt / self.domAltMaxSpeed
            domAzSlewTime  = domDeltaAz  / self.domAzMaxSpeed
            #print "domAlt1", domAltSlewTime
            #print "domAz1", domAzSlewTime

            totDomTime1 = max(domAltSlewTime, domAzSlewTime)

            domDeltaAlt = deltaAlt - 2 * (self.domSlitDiam/2 - self.fovWidth/2)
            domDeltaAz  = deltaAz
            domAltSlewTime = domDeltaAlt / self.domAltMaxSpeed
            domAzSlewTime  = domDeltaAz  / self.domAzMaxSpeed
            #print "domAlt2", domAltSlewTime
            #print "domAz2", domAzSlewTime
            totDomTime2 = max(domAltSlewTime, domAzSlewTime)

            totDomTime = min(totDomTime1, totDomTime2)

        #print "slew alt/az", altSlewTime, azSlewTime
        #print
        #print "tel Alt", telAltSlewTime
        #print "tel Az", telAzSlewTime
        totTelTime = max(telAltSlewTime, telAzSlewTime) + self.settleTime
        slewTime = max(totTelTime, totDomTime)

        # include filter change time if necessary
        if filter1 != filter2:
            slewTime = max(slewTime, self.filterChangeTime)

        return slewTime
