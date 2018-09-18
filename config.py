from __future__ import division

import numpy as np
import time
from lsst.sims.speedObservatory import Telescope

# number of nights to run the simulation for
surveyNumNights = 365 * 10

# start time of the survey as a unix timestamp
# this is january 1 2022 at midnight
surveyStartTime = 1640995200

# minimum and maximum declinations to observe at
# -90 to 13 gets about equal area in the North and Southeast
maxDec = np.radians(13)
minDec = np.radians(-90)

# -75 to 15 is minion_1020 (pan-starrs like)
#maxDec = np.radians(15)
#minDec = np.radians(-75)

# -62 to 2.5 is minion_1012 (wfd only)
#maxDec = np.radians(3)
#minDec = np.radians(-64)

# Approximate width in RA of scans in the zenith dec band
# 30 degrees seems reasonable
EScanWidth = np.radians(30)

# buffer zone near the zenith to the N/S/E 
# TODO these probably shouldn't change, so maybe shouldn't be config vars
# otoh zenithBufferOffset is a free parameters that'll change the survey
# properties
# TODO these names are terrible I have no idea what they are w/o the comment
tel = Telescope()
# zenithBuffer is half of the width (in dec) of the zenith dec band
zenithBuffer = 2 * tel.fovWidth
# zenithBufferOffset is the distance the zenith dec band is offset
# in RA so that you avoid the zenith avoidance zone
zenithBufferOffset = 4 * tel.fovWidth

# density of the tiling in units of pointings / steradian
# this may be the approximate value depending on whether the Tiling
# subclass can achieve arbitrary densities or not (e.g. ThomsonTiling can't)
tilingDensity = (1 / 12) * (180/np.pi)**2

# don't observe when the cloud is more than maxCloudCover cloudy
maxCloudCover = 0.7

# don't observe in u if the moon is above moonUMaxAlt and if the moon's phase
# is above moonUMaxPhase
moonUMaxAlt = np.radians(0)
moonUMaxPhase = 0.3

# similar for g
moonGMaxAlt = np.radians(20)
moonGMaxPhase = 0.6

# whether or not to include deep drilling fields
useDD = True

# number of seconds of exposure time for DD visits
DDExpTime = 60 * 60

# total number of seconds of exposure time for WFD visits
# so 2 15 second exposures would correspond to 30 here
WFDExpTime = 30

# exposure time in twilight
xTwilExpTime = 5

# number of exposures that are taken per visit
numExposures = 2
xTwilNumExposures = 1

# number of seconds of overhead needed per visit
# (shutter and potentially intermediate readout times)
visitOverheadTime = (numExposures - 1) * 3 + 1
xTwilVisitOverheadTime = (xTwilNumExposures - 1) * 3 + 1
visitOverheadTime = xTwilVisitOverheadTime

# whether to use the "relaxed" dome model in the sims.speedObservatory
# telescope model
laxDome = False

# constants (the values don't matter)
NORTH = 0
SOUTH = 1
EAST = 2
WEST = 3
SOUTHEAST = 4
