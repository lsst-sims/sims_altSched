from __future__ import division

import numpy as np
import time
from Telescope import Telescope

surveyNumNights = 365 * 10
# this is january 1 2022 at midnight
surveyStartTime = 1640995200

maxDec = np.radians(30)
minDec = np.radians(-90)

# buffer zone near the zenith to the N/S/E 
# TODO these probably shouldn't change, so probably shouldn't be config vars
tel = Telescope()
zenithBuffer = 2 * tel.fovWidth
zenithBufferOffset = 4 * tel.fovWidth
