from __future__ import division

import numpy as np
import time
from Telescope import Telescope

surveyNumNights = 365 * 10
# this is january 1 2022 at midnight
surveyStartTime = 1640995200

# 6.2 degrees is where the area in the south and east
# equals the area in the north
maxDec = np.radians(6.2)
minDec = np.radians(-90)

# buffer zone near the zenith to the N/S/E 
# TODO these probably shouldn't change, so probably shouldn't be config vars
tel = Telescope()
zenithBuffer = 2 * tel.fovWidth
zenithBufferOffset = 4 * tel.fovWidth
