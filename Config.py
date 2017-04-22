from __future__ import division

import numpy as np
import time
import Telescope

surveyNumNights = 365 * 10
# this is january 1 2022 at midnight
surveyStartTime = 1640995200

maxDec = np.radians(30)
minDec = np.radians(-90)

# buffer zone near the zenith to the N/S/E 
# TODO these probably shouldn't change, so probably shouldn't be config vars
zenithBuffer = 2 * Telescope.fovWidth
zenithBufferOffset = 4 * Telescope.fovWidth
