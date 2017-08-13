from __future__ import division

import numpy as np
import time
from lsst.sims.speedObservatory import Telescope

surveyNumNights = 365 * 10
# this is january 1 2022 at midnight
surveyStartTime = 1640995200

# 6.2 degrees should be where the area in the south and east
# equals the area in the north but apparently not?
maxDec = np.radians(11)
minDec = np.radians(-90)

# Approximate width in RA of easterly scans
# 30 degrees seems reasonable
EScanWidth = np.radians(30)

# buffer zone near the zenith to the N/S/E 
# TODO these probably shouldn't change, so maybe shouldn't be config vars
# otoh zenithBufferOffset is a free parameters that'll change the survey
# properties
# TODO these names are terrible I have no idea what they are w/o the comment
tel = Telescope()
# zenithBuffer is half of the width of the zenith dec band
zenithBuffer = 2 * tel.fovWidth
# zenithBufferOffset is the distance the zenith dec band is offset
# in RA so that you avoid the zenith avoidance zone
zenithBufferOffset = 4 * tel.fovWidth

# don't observe when the cloud is more than maxCloudCover cloudy
maxCloudCover = 0.7

moonUMaxAlt = np.radians(0)
moonUMaxPhase = 0.3

moonGMaxAlt = np.radians(20)
moonGMaxPhase = 0.6

# number of seconds of overhead needed per visit
# (shutter and potentially intermediate readout times)
visitOverheadTime = 4

NORTH = 0
SOUTH = 1
EAST = 2
WEST = 3
SOUTHEAST = 4
