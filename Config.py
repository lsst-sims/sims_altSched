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
# TODO these probably shouldn't change, so probably shouldn't be config vars
tel = Telescope()
zenithBuffer = 2 * tel.fovWidth
zenithBufferOffset = 4 * tel.fovWidth

# don't observe when the cloud is more than maxCloudCover cloudy
maxCloudCover = 0.7

# number of seconds of overhead needed per visit
# (shutter and potentially intermediate readout times)
visitOverheadTime = 4
