from __future__ import print_function, division

import numpy as np
import sys
import csv
import healpy as hp
from datetime import datetime
import palpy
import time
import pandas
import sqlite3

from lsst.sims.speedObservatory import sky, utils
from lsst.sims.speedObservatory import Telescope
#from lsst.sims.skybrightness import SkyModel
from lsst.sims.skybrightness_pre import SkyModelPre
from lsst.sims.utils import raDec2Hpid
from lsst.sims.utils import m5_flat_sed
from lsst.sims.ocs.kernel.time_handler import TimeHandler
from lsst.sims.ocs.environment.cloud_model import CloudModel
from lsst.sims.ocs.environment.seeing_model import SeeingModel
from lsst.sims.ocs.configuration.instrument import Filters
from lsst.sims.ocs.configuration import Environment
import config

inFilename = sys.argv[1]
outFilename = sys.argv[2]

# input is a CSV file with
# mjd, prop, ra, dec, filter

# output is as close as we can get to OpSim's SummaryAllProps table:
# observationId (int)
# night (int)
# observationStartTime (float, UTC, units=seconds)
# observationStartMJD (float, units=days)
# observationStartLST (float, units=degrees)
# filter (string(1))
# proposalId (int)
# fieldRA (units=degrees)
# fieldDec (units=degrees)
# angle (position angle, units=degrees)
# altitude (units=degrees)
# azimuth (units=degrees)
# numExposures (int)
# visitTime (float, exp+overhead, units=seconds)
# visitExposureTime (float, exp, units=seconds)
# airmass (float)
# skyBrightness (float)
# cloud (float)
# seeingFwhm500 (float)
# seeingFwhmGeom (float)
# seeingFwhmEff (float)
# fiveSigmaDepth (float, units=magnitudes)
# moonRA (float, units=degrees)
# moonDec (float, units=degrees)
# moonAlt (float, units=degrees)
# moonAz (float, units=degrees)
# moonDistance (float, units=degrees)
# moonPhase (float)
# sunRA (float, units=degrees)
# sunDec (float, units=degrees)
# sunAlt (float, units=degrees)
# sunAz (float, units=degrees)
# solarElong (float, units=degrees)
# slewTime (float, units=seconds)
# slewDistance (float, units=degrees)
# paraAngle (float, units=degrees)
# rotTelPos (float, units=degrees)
# rotSkyPos (float, units=degrees)

########################################

# first read in the input
inData = []
with open(inFilename, "r") as inFile:
    reader = csv.DictReader(inFile)
    for row in reader:
        inData.append((float(row["mjd"]),
                       row["prop"],
                       float(row["ra"]),
                       float(row["dec"]),
                       row["filter"]
                      ))
names = ["mjd", "prop",       "ra",  "dec", "filter"]
types = [float, (np.str_, 3), float, float, (np.str_, 1)]
inDtype = list(zip(names, types))
inData = np.array(inData, dtype=np.dtype(inDtype))
nRows = inData.shape[0]

inputType = "opsim"
assert(inputType == "alt-sched" or inputType == "opsim")
print("Using input type {}".format(inputType))

if inputType == "alt-sched":
    # dd visits are treated as single long exposures
    props = ["wfd", "dd"]
    propVisitTimes = [config.WFDExpTime + config.visitOverheadTime,
                      config.DDExpTime]
    propVisitExpTimes = [config.WFDExpTime, config.DDExpTime]
    propNames = ["WideFastDeep", "Deep Drilling"]

elif inputType == "opsim":
    # all types of visits in opsim are 30 + 4 seconds
    props = ["wfd", "dd", "gal", "ecl", "scp"]
    propVisitTimes = [34] * 5
    propVisitExpTimes = [30] * 5
    propNames = ["WideFastDeep", "Deep Drilling", "Galactic Plane",
                 "North Ecliptic Spur", "South Celestial Pole"]

propIds = np.arange(len(props))

outDtype = [
            ("observationId", int),
            ("night", int),
            ("observationStartTime", float),
            ("observationStartMJD", float),
            ("observationStartLST", float),
            ("filter", (np.str_, 1)),
            ("proposalId", (int)),
            ("fieldRA", float),
            ("fieldDec", float),
            ("angle", float),
            ("altitude", float),
            ("azimuth", float),
            ("numExposures", int),
            ("visitTime", float),
            ("visitExposureTime", float),
            ("airmass", float),
            ("skyBrightness", float),
            ("cloud", float),
            ("seeingFwhm500", float),
            ("seeingFwhmGeom", float),
            ("seeingFwhmEff", float),
            ("fiveSigmaDepth", float),
            ("moonRA", float),
            ("moonDec", float),
            ("moonAlt", float),
            ("moonAz", float),
            ("moonDistance", float),
            ("moonPhase", float),
            ("sunRA", float),
            ("sunDec", float),
            ("sunAlt", float),
            ("sunAz", float),
            ("solarElong", float),
            ("slewTime", float),
            ("slewDistance", float),
            ("paraAngle", float),
            ("rotTelPos", float),
            ("rotSkyPos", float)
           ]

outData = np.zeros(nRows, dtype=np.dtype(outDtype))

class Timer:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        interval = time.clock() - self.start
        print("{:<21}: {:2.2f} secs".format(self.name, interval))

### observationId ###
with Timer("observationId") as t:
    outData["observationId"] = np.arange(nRows)

"""
### night ###
# this is a stupid/slow way to calculate night. Instead, use
# the less stupid way below
with Timer("night") as t:
    # it's expensive to call sky.twilStart/End, so we only calculate
    # it config.surveyNumNights / nightStride times and then estimate
    # the twilStart/End times for the remaining nights

    nightRange = np.arange(config.surveyNumNights + 1)
    twilStarts = []
    twilEnds = []

    nightNum = 0
    nightStride = 20
    # 60 minute buffer to subtract off twilStarts and to add to twilEnds
    # we just use these start/end times to calculate night numbers for each
    # observation, so it doesn't matter if we pretend the night is longer than
    # it is but it does matter if we mistakenly assume it's shorter than it is
    timeBuffer = 3600
    while nightNum < config.surveyNumNights + 1:
        # get night start/end times for nightNum
        twilStart = sky.twilStart(config.surveyStartTime, nightNum)
        twilEnd = sky.twilEnd(config.surveyStartTime, nightNum)

        # save start/end times
        twilStarts.append(twilStart - timeBuffer)
        twilEnds.append(twilEnd + timeBuffer)

        # calculate the time delta between the start/end for nightNum + 1
        # and the start/end for nightNum
        nextTwilStart = sky.twilStart(config.surveyStartTime, nightNum + 1)
        startInterval = nextTwilStart - twilStart
        nextTwilEnd = sky.twilEnd(config.surveyStartTime, nightNum + 1)
        endInterval = nextTwilEnd - twilEnd

        # use this time delta to estimate the following nightStride nights'
        # start/end times
        for i in range(1, nightStride):
            twilStarts.append(twilStart + i * startInterval - timeBuffer)
            twilEnds.append(twilEnd + i * startInterval + timeBuffer)
        nightNum += nightStride

    twilStarts = np.array(twilStarts)
    twilEnds = np.array(twilEnds)
    twilStarts = utils.unix2mjd(twilStarts)
    twilEnds = utils.unix2mjd(twilEnds)
    nights = np.searchsorted(twilStarts, inData["mjd"], side="right") - 1
    # check that the inData mjds are between their nights' start and end times
    print("inData < starts", (twilStarts[nights] > inData["mjd"]).sum())
    print("inData > ends", (twilEnds[nights] < inData["mjd"]).sum())

    bad = np.where((twilStarts[nights] > inData["mjd"]) | (inData["mjd"] > twilEnds[nights]))[0]
    print((inData["mjd"][bad] - twilEnds[nights[bad]]) * 3600*24)
    assert(((twilStarts[nights] <= inData["mjd"]) &
            (inData["mjd"] <= twilEnds[nights])).all())
    outData["night"] = nights
"""

### night ###
with Timer("night_check") as t:
    timeBuffer = 1800
    if inputType == "alt-sched":
        startTime = config.surveyStartTime
    elif inputType == "opsim":
        startTime = 1664500000
        #startTime = 1664582400 # for astro-lsst
        #startTime = config.surveyStartTime - 3600*24
    nightRange = np.arange(config.surveyNumNights + 3)
    twilStarts = np.array([sky.twilStart(startTime, night)
                           for night in nightRange]) - timeBuffer
    twilEnds   = np.array([sky.twilEnd(startTime, night)
                           for night in nightRange]) + timeBuffer
    twilStarts = utils.unix2mjd(twilStarts)
    twilEnds = utils.unix2mjd(twilEnds)
    nights = np.searchsorted(twilStarts, inData["mjd"], side="right") - 1
    # check that the inData mjds are between their nights' start and end times
    try:
        assert(((twilStarts[nights] <= inData["mjd"]) &
                (inData["mjd"] <= twilEnds[nights])).all())
    except AssertionError as e:
        print(twilStarts[nights])
        print(inData["mjd"])
        print(twilEnds[nights])
        print(nights)
        print("twilStarts > inData['mjd']", (twilStarts[nights] > inData["mjd"]).sum())
        print("twilEnds < inData['mjd']", (twilEnds[nights] < inData["mjd"]).sum())
        overRun = np.where(twilEnds[nights] < inData["mjd"])
        print(np.unique(nights[overRun]))
        print(twilStarts[nights[overRun]])
        print(inData["mjd"][overRun])
        print(twilEnds[nights[overRun]])
        raise e

    outData["night"] = nights



### observationStartTime ###
with Timer("observationStartTime") as t:
    outData["observationStartTime"] = utils.mjd2unix(inData["mjd"])


### observationStartMJD ###
with Timer("observationStartMJD") as t:
    outData["observationStartMJD"] = inData["mjd"]


### observationStartLST ###
with Timer("observationStartLST") as t:
    telLongitude = Telescope.longitude
    lstRad = sky.unix2lst(telLongitude, utils.mjd2unix(inData["mjd"]))
    outData["observationStartLST"] = np.degrees(lstRad)


### filter ###
with Timer("filter") as t:
    outData["filter"] = inData["filter"]

### proposalId ###
with Timer("proposalId") as t:
    # put an invalid value to start
    outData["proposalId"] = min(propIds) - 1

    for propId, prop in zip(propIds, props):
        outData["proposalId"][np.where(inData["prop"] == prop)] = propId

    # make sure no invalid values remain
    assert(min(outData["proposalId"]) >= min(propIds))

### ra ###
with Timer("fieldRA") as t:
    outData["fieldRA"] = np.degrees(inData["ra"])


### dec ###
with Timer("fieldDec") as t:
    outData["fieldDec"] = np.degrees(inData["dec"])


### angle ###
with Timer("angle") as t:
    outData["angle"] = np.zeros(nRows)


### altitude ###
with Timer("altitude") as t:
    alt, az = sky.radec2altaz(inData["ra"], inData["dec"], utils.mjd2unix(inData["mjd"]))
    outData["altitude"] = np.degrees(alt)


### azimuth ###
with Timer("azimuth") as t:
    outData["azimuth"] = np.degrees(az)


### numExposures ###
with Timer("numExposures") as t:
    outData["numExposures"] = config.numExposures


### visitTime ###
### visitExposureTime ###
with Timer("visit(Exp)Time") as t:
    # keep track of how many props are assigned to each visit (just for assertion)
    assigned = np.zeros(nRows)

    for propId, prop in zip(propIds, props):
        thisProp = inData["prop"] == prop
        assigned[thisProp] += 1
        outData["visitTime"][thisProp] = propVisitTimes[propId]
        outData["visitExposureTime"][thisProp] = propVisitExpTimes[propId]

    # make sure all visits get exactly one prop
    assert(assigned.min() == assigned.max() == 1)


### airmass ###
with Timer("airmass") as t:
    outData["airmass"] = 1 / np.sin(alt)



### skyBrightness ###
# this takes about 1000 secs (with files on SSD)
with Timer("skyBrightness") as t:
    sm = SkyModelPre(data_path="/home/doctor/tmp/healpix", verbose=False)
    full = sm.returnMags(inData["mjd"][0], airmass_mask=False, planet_mask=False,
                                           moon_mask=False, zenith_mask=False)
    nside = hp.npix2nside(full["r"].size)
    skyBrightnesses = np.zeros(nRows)
    for i, obs in enumerate(inData):
        indx = raDec2Hpid(nside, np.degrees(obs["ra"]), np.degrees(obs['dec']))
        skyBrightnesses[i] = \
                sm.returnMags(obs['mjd'], indx=[indx], airmass_mask=False,
                              planet_mask=False, moon_mask=False,
                              zenith_mask=False, extrapolate=True)[obs['filter']]
        if outData["night"][i] == 197 and outData["filter"][i] == "u":
            print(obs["mjd"], obs["ra"], obs["dec"], indx, skyBrightnesses[i])

        print("Progress:", "{:3.4f}%".format(i / nRows * 100), end="\r")
    del sm
    outData["skyBrightness"] = skyBrightnesses


"""
# this takes about 3600 secs
with Timer("skyBrightness") as t:
    skyModel = SkyModel(mags=True)
    for i, inDatum in enumerate(inData):
        skyModel.setRaDecAltAzMjd(np.array([inDatum["ra"]]),
                                  np.array([inDatum["dec"]]),
                                  np.array([np.radians(outData["altitude"][i])]),
                                  np.array([np.radians(outData["azimuth"][i])]),
                                  inDatum["mjd"])
        
        outData["skyBrightness"][i] = skyModel.returnMags()[inDatum["filter"]]
        print("Progress:", "{:3.4f}%   ".format(i / nRows * 100), end="\r")
"""

### cloud ###
deltaTs = utils.mjd2unix(inData["mjd"]) - config.surveyStartTime
with Timer("cloud") as t:
    # make a TimeHandler for the cloud database
    dateFormat = "%Y-%m-%d"
    startDatetime = datetime.utcfromtimestamp(config.surveyStartTime)
    timeHandler = TimeHandler(datetime.strftime(startDatetime, dateFormat))
    cloudModel = CloudModel(timeHandler)
    # load the cloud database
    cloudModel.initialize()
    cloudCover = [cloudModel.get_cloud(deltaT) for deltaT in deltaTs]
    outData["cloud"] = cloudCover


### seeingFwhm500 ###
with Timer("seeingFwhm500") as t:
    seeingModel = SeeingModel(timeHandler)
    environment = Environment()
    filters = Filters()
    seeingModel.initialize(environment, filters)
    fwhm500s = np.zeros(nRows)
    fwhmGeoms = np.zeros(nRows)
    fwhmEffs = np.zeros(nRows)
    for i, deltaT, f, airmass in zip(np.arange(nRows), deltaTs,
                                     inData["filter"], outData["airmass"]):
        fwhm500, fwhmGeom, fwhmEff = seeingModel.calculate_seeing(deltaT, f, airmass)
        fwhm500s[i] = fwhm500
        fwhmGeoms[i] = fwhmGeom
        fwhmEffs[i] = fwhmEff

    outData["seeingFwhm500"] = fwhm500s


### seeingFwhmGeom ###
with Timer("seeingFwhmGeom") as t:
    outData["seeingFwhmGeom"] = fwhmGeoms


### seeingFwhmEff ###
with Timer("seeingFwhmEff") as t:
    outData["seeingFwhmEff"] = fwhmEffs


### fiveSigmaDepth ###
with Timer("fiveSigmaDepth") as t:
    for f in Telescope.filters:
        relevantRows = np.where(inData["filter"] == f)
        fiveSigmaDepth = m5_flat_sed(f,
                                     outData["skyBrightness"][relevantRows],
                                     outData["seeingFwhmEff"][relevantRows],
                                     outData["visitExposureTime"][relevantRows],
                                     outData["airmass"][relevantRows])
        outData["fiveSigmaDepth"][relevantRows] = fiveSigmaDepth


### moonRa ###
moonRas = np.zeros(nRows)
moonDecs = np.zeros(nRows)
with Timer("moonRa") as t:
    for i, unix in enumerate(utils.mjd2unix(inData["mjd"])):
        moonRa, moonDec = sky.radecOfMoon(unix)
        moonRas[i] = moonRa
        moonDecs[i] = moonDec
    outData["moonRA"] = np.degrees(moonRas)


### moonDec ###
with Timer("moonDec") as t:
    outData["moonDec"] = np.degrees(moonDecs)


### moonAlt ###
moonAlt = None
moonAz = None
with Timer("moonAlt") as t:
    moonAlt, moonAz = sky.radec2altaz(moonRas, moonDecs, utils.mjd2unix(inData["mjd"]))
    outData["moonAlt"] = np.degrees(moonAlt)


### moonAz ###
with Timer("moonAz") as t:
    outData["moonAz"] = np.degrees(moonAz)


### moonDistance ###
with Timer("moonDistance") as t:
    moonDist = palpy.dsepVector(inData["ra"], inData["dec"], moonRas, moonDecs)
    outData["moonDistance"] = np.degrees(moonDist)

### moonPhase ###
with Timer("moonPhase") as t:
    outData["moonPhase"] = [sky.phaseOfMoon(unix)
                            for unix in utils.mjd2unix(inData["mjd"])]


### sunRA ###
sunRas = np.zeros(nRows)
sunDecs = np.zeros(nRows)
with Timer("sunRA") as t:
    for i, unix in enumerate(utils.mjd2unix(inData["mjd"])):
        sunRa, sunDec = sky.radecOfSun(unix)
        sunRas[i] = sunRa
        sunDecs[i] = sunDec
    outData["sunRA"] = np.degrees(sunRas)

### sunDec ###
with Timer("sunDec") as t:
    outData["sunDec"] = np.degrees(sunDecs)


### sunAlt ###
sunAlt = None
sunAz = None
with Timer("sunAlt") as t:
    sunAlt, sunAz = sky.radec2altaz(sunRas, sunDecs, utils.mjd2unix(inData["mjd"]))
    outData["sunAlt"] = np.degrees(sunAlt)


### sunAz ###
with Timer("sunAz") as t:
    outData["sunAz"] = np.degrees(sunAz)


### solarElong ###
with Timer("solarElong") as t:
    sunDist = palpy.dsepVector(inData["ra"], inData["dec"], sunRas, sunDecs)
    outData["solarElong"] = np.degrees(sunDist)


### slewTime ###
with Timer("slewTime") as t:
    # make a default sims_speedObservatory telescope

    tel = Telescope()
    prevMjd = 0
    prevAlt = 0
    prevAz = 0
    prevF = "u"

    alts = np.radians(outData["altitude"])
    azes = np.radians(outData["azimuth"])
    for i, mjd, alt, az, f in zip(np.arange(nRows), inData["mjd"],
                                  alts, azes, inData["filter"]):
        if mjd - prevMjd > 1 / 24:
            slewTime = 0
        else:
            slewTime = tel.calcSlewTime(prevAlt, prevAz, prevF,
                                        alt, az, f,
                                        laxDome = config.laxDome)
        outData["slewTime"][i] = slewTime
        prevMjd = mjd
        prevAlt = alt
        prevAz = az
        prevF = f


### slewDistance ###
with Timer("slewDistance") as t:
    alt = np.radians(outData["altitude"])
    az = np.radians(outData["azimuth"])
    slewDist = palpy.dsepVector(az[1:], alt[1:],
                                az[:-1], alt[:-1])
    outData["slewDistance"][1:] = np.degrees(slewDist)
    outData["slewDistance"][outData["slewTime"]==0] = 0

# paraAngle (float, units=degrees)
# rotTelPos (float, units=degrees)
# rotSkyPos (float, units=degrees)

#######################################3
exit()

# now output to database
conn = sqlite3.connect(outFilename)

dataFrame = pandas.DataFrame(outData)
dataFrame.to_sql("SummaryAllProps", conn)

# make Proposal table:
propDtype = [("propId", int), ("propName", (np.str_, 32))]
propData = np.zeros(len(propIds), dtype=np.dtype(propDtype))

propData["propId"] = propIds
propData["propName"] = propNames
propDataFrame = pandas.DataFrame(propData)
propDataFrame.to_sql("Proposal", conn)
