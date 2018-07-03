#!/usr/bin/env python

from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
import os
import argparse
import copy
import numpy as np
import warnings
import matplotlib
# Set matplotlib backend (to create plots where DISPLAY is not set).
matplotlib.use('Agg')
import matplotlib.cm as cm

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.plots as plots
import lsst.sims.maf.utils as utils

class SurveyTimeStacker(stackers.BaseStacker):
    def __init__(self, visitTimeCol="visitTime", slewTimeCol="slewTime"):
        self.visitTimeCol = visitTimeCol
        self.slewTimeCol = slewTimeCol
        self.colsReq = [visitTimeCol, slewTimeCol]
        self.colsAdded = ["surveyTime"]
        self.units = ["seconds"]

    def _run(self, simData):
        simData["surveyTime"] = simData[self.visitTimeCol] + simData[self.slewTimeCol]
        return simData

class NormSumMetric(metrics.BaseMetric):
    def __init__(self, col, normValue, metricName="Normalized Sum", **kwargs):
        super(NormSumMetric, self).__init__(col=col, metricName=metricName, **kwargs)
        self.col = col
        self.normValue = normValue
    def run(self, dataSlice, slicePoint=None):
        return np.sum(dataSlice[self.col]) / self.normValue



def makeBundleList(dbFile, runName=None, benchmark='design'):

    seeingCol = 'seeingFwhmEff'
    benchmarkSeeing = 'FWHMeff'
    latCol = "fieldDec"
    lonCol = "fieldRA"
    nside=64

    # List to hold everything we're going to make
    bundleList = []

    # Connect to the databse
    opsimdb = db.OpsimDatabase(dbFile)
    propCol = 'proposalId'

    if runName is None:
        runName = os.path.basename(dbFile).replace('_sqlite.db', '')
        runName = runName.replace('.db', '')

    # Fetch the proposal ID values from the database
    propids, propTags = opsimdb.fetchPropInfo()
    DDpropid = propTags['DD']
    WFDpropid = propTags['WFD']

    # Fetch the telescope location from config
    lat, lon, height = opsimdb.fetchLatLonHeight()

    # Construct a WFD SQL where clause so multiple propIDs can query by WFD:
    wfdWhere = opsdb.createSQLWhere('WFD', propTags)
    print('#FYI: WFD "where" clause: %s' % (wfdWhere))
    ddWhere = opsdb.createSQLWhere('DD', propTags)
    print('#FYI: DD "where" clause: %s' % (ddWhere))

    # Set up benchmark values, scaled to length of opsim run. These are applied to 'all' and 'WFD' plots.
    runLength = opsimdb.fetchRunLength()
    if benchmark == 'requested':
        # Fetch design values for seeing/skybrightness/single visit depth.
        benchmarkVals = utils.scaleBenchmarks(runLength, benchmark='design')
        # Update nvisits with requested visits from config files.
        benchmarkVals['nvisits'] = opsimdb.fetchRequestedNvisits(propId=WFDpropid)
        # Calculate expected coadded depth.
        benchmarkVals['coaddedDepth'] = utils.calcCoaddedDepth(
            benchmarkVals['nvisits'], benchmarkVals['singleVisitDepth'])
    elif (benchmark == 'stretch') or (benchmark == 'design'):
        # Calculate benchmarks for stretch or design.
        benchmarkVals = utils.scaleBenchmarks(runLength, benchmark=benchmark)
        benchmarkVals['coaddedDepth'] = utils.calcCoaddedDepth(
            benchmarkVals['nvisits'], benchmarkVals['singleVisitDepth'])
    else:
        raise ValueError(
            'Could not recognize benchmark value %s, use design, stretch or requested.' % (benchmark))
    # Check that nvisits is not set to zero (for very short run length).
    for f in benchmarkVals['nvisits']:
        if benchmarkVals['nvisits'][f] == 0:
            print('Updating benchmark nvisits value in %s to be nonzero' % (f))
            benchmarkVals['nvisits'][f] = 1


    # Set values for min/max range of nvisits for All/WFD and DD plots. These are somewhat arbitrary.
    nvisitsRange = {}
    nvisitsRange['all'] = {'u': [30, 110], 'g': [50, 160], 'r': [100, 310],
                           'i': [90, 270], 'z': [100, 330], 'y': [50, 300]}

    # Scale these nvisit ranges for the runLength.
    scale = runLength / 10.0
    for prop in nvisitsRange:
        for f in nvisitsRange[prop]:
            for i in [0, 1]:
                nvisitsRange[prop][f][i] = int(np.floor(nvisitsRange[prop][f][i] * scale))

    # Filter list, and map of colors (for plots) to filters.
    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    colors = {'u': 'm', 'g': 'b', 'r': 'g', 'i': 'y', 'z': 'r', 'y': 'k'}
    filtorder = {'u': 1, 'g': 2, 'r': 3, 'i': 4, 'z': 5, 'y': 6}

    slicermetadata = ''

    ###
    # Configure some standard summary statistics dictionaries to apply to appropriate metrics.
    standardStats = [metrics.MeanMetric(),
                     metrics.RmsMetric(), metrics.MedianMetric(), metrics.CountMetric(),
                     metrics.NoutliersNsigmaMetric(metricName='N(+3Sigma)', nSigma=3),
                     metrics.NoutliersNsigmaMetric(metricName='N(-3Sigma)', nSigma=-3.)]

    rangeStats = [metrics.PercentileMetric(metricName='25th%ile', percentile=25),
                  metrics.PercentileMetric(metricName='75th%ile', percentile=75),
                  metrics.MinMetric(), metrics.MaxMetric()]

    allStats = copy.deepcopy(standardStats)
    allStats.extend(rangeStats)

    # Standardize a couple of labels (for ordering purposes in showMaf).
    summaryGroup = 'A: Summary'
    nvisitGroup = 'B: NVisits'
    nvisitPerPropGroup = 'C: NVisits (per prop)'
    surveyTimeGroup = 'D: Survey Time'
    hourangleGroup = 'E: Hour Angle'
    distanceGroup = 'F: Distance to Sun/Moon'
    hourglassGroup = 'G: Hourglass'
    filterGroup = 'H: Filter Changes'
    slewGroup = 'I: Slew'

    # Fetch the total number of visits (to create fraction for number of visits per proposal)
    totalNVisits = opsimdb.fetchNVisits()
    totalSlewN = opsimdb.fetchTotalSlewN()

    # Set up an object to hold all the bundles that will be merged together
    healpixHistPlot = plots.HealpixHistogram()
    mergedHistDict = {}
    healPlots = ['nvisits', 'fullRangeHA', 'meanHourAngle']
    oneDPlots = ['moonDistance', 'hourAngle'] 

    # set plotters for merged histograms (healpix hist or oned hist)
    for plotName in healPlots:
        mergedHistDict[plotName] = plots.PlotBundle(plotFunc=healpixHistPlot)
    for plotName in oneDPlots:
        mergedHistDict[plotName] = plots.PlotBundle(plotFunc=plots.OneDBinnedData())

    # Metrics calculating values across the sky.
    # Loop over a set of standard analysis metrics, for All Proposals, WFD only, and DD only.

    prop = "WFD"
    for f in filters:
        propCaption = 'for all WFD proposals'
        metadata = '%s band, WFD' % (f) + slicermetadata
        sqlconstraint = 'filter = "%s" and %s' % (f, wfdWhere)
        nvisitsMin = nvisitsRange['all'][f][0]
        nvisitsMax = nvisitsRange['all'][f][1]
        mag_zp = benchmarkVals['coaddedDepth'][f]

        # Make a new slicer for each metric since they can get setup with different fields later
        slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)

        # Configure the metrics to run for this sql constraint (all proposals/wfd and filter combo).

        # Count the total number of visits.
        metric = metrics.CountMetric(col='observationStartMJD', metricName='NVisits')
        plotDict = {'xlabel': 'Number of Visits',
                    'xMin':     nvisitsMin, 'xMax':     nvisitsMax,
                    'colorMin': nvisitsMin, 'colorMax': nvisitsMax,
                    'binsize': 5}
        summaryStats = allStats
        displayDict = {'group': nvisitGroup, 'order': filtorder[f], "subgroup": "per filter",
                       'caption': 'Number of visits in filter %s, %s.' % (f, propCaption)}
        histMerge = {'color': colors[f], 'label': '%s' % (f),
                     'xMin': nvisitsMin, 'xMax': nvisitsMax,
                     'binsize': 5, 'legendloc': 'upper right'}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata,
                                            summaryMetrics=summaryStats)
        mergedHistDict['nvisits'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Calculate the mean of the hour angle.
        metric = metrics.MeanMetric(col='HA')
        plotDict = {'xMin':     -5, 'xMax':     5,
                    'colorMin': -5, 'colorMax': 5,
                    'binsize': 0.05}
        displayDict = {'group': hourangleGroup, 'order': filtorder[f], "subgroup": "mean",
                       'caption': 'Mean of the Hour Angle in filter %s, %s.'
                       % (f, propCaption)}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata,
                                            summaryMetrics=summaryStats)
        histMerge = {'color': colors[f], 'label': '%s' % (f),
                     'xlabel': 'Mean Hour Angle (Hours)',
                     'xMin': -5, 'xMax': 5,
                     'binsize': .05, 'legendloc': 'upper right'}
        mergedHistDict['meanHourAngle'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Calculate the Full Range of the hour angle.
        metric = metrics.FullRangeMetric(col='HA')
        plotDict = {'xMin':     0, 'xMax':     12,
                    'colorMin': 0, 'colorMax': 12,
                    'binsize': 0.1}
        displayDict = {'group': hourangleGroup, 'order': filtorder[f], "subgroup": "full range",
                       'caption': 'Full Range of the Hour Angle in filter %s, %s.'
                       % (f, propCaption)}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata,
                                            summaryMetrics=summaryStats)
        histMerge = {'color': colors[f], 'label': '%s' % (f),
                     'xlabel': 'Full Hour Angle Range',
                     'binsize': 0.1,
                     'legendloc': 'upper right'}
        mergedHistDict['fullRangeHA'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)


    # number of visits in all filters together, WFD only.
    sqlconstraint = wfdWhere
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
    metric = metrics.CountMetric(col='observationStartMJD', metricName='NVisits')
    metadata = 'All filters, WFD'
    plotDict = {'xlabel': 'Number of Visits', 'binsize': 5, 'cumulative': False,
                'xMin': 500, 'xMax': 1500,
                "colorMin": 500, "colorMax": 1500}
    summaryStats = allStats
    displayDict = {'group': nvisitPerPropGroup, 'subgroup': 'All Filters', 'order': 0,
                   'caption': 'Number of visits all filters, WFD only'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    # Number of visits in all filters together, DD only.
    sqlconstraint = ddWhere
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
    metric = metrics.CountMetric(col='observationStartMJD', metricName='NVisits')
    metadata = 'All filters, DD'
    plotDict = {'xlabel': 'Number of Visits', 'binsize': 5, 'cumulative': False,
                'xMin': 0, 'xMax': 1000,
                "colorMin": 0, "colorMax": 1000}
    summaryStats = allStats
    displayDict = {'group': nvisitPerPropGroup, 'subgroup': 'All Filters', 'order': 1,
                   'caption': 'Number of visits all filters, DD only'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    sqlconstraint = wfdWhere
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)

    stackerList = [SurveyTimeStacker()]
    metric = NormSumMetric(col="surveyTime", normValue=60, metricName="Survey Time (mins)")
    metadata = "All filters, WFD"
    plotDict = {"xlabel": "Survey Time (minutes)",
                "xMin": 500, "xMax": 800, "binsize": 10,
                "colorMin": 500, "colorMax": 800}
    summaryStats = allStats
    displayDict = {"group": surveyTimeGroup, "subgroup": "WFD",
                   "caption": "Amount of survey time (visit time plus slew time)" + \
                              " spent on each sky pixel"}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName,
                                        metadata=metadata, summaryMetrics=summaryStats,
                                        stackerList=stackerList)
    bundleList.append(bundle)

    # End of all-sky metrics.

    # Hourglass metric.
    hourSlicer = slicers.HourglassSlicer()
    # Calculate Filter Hourglass plots per year (split to make labelling easier).
    yearDates = list(range(0, int(round(365 * runLength)) + 365, 365))
    for i in range(len(yearDates) - 1):
        sqlconstraint = 'night > %i and night <= %i' % (yearDates[i], yearDates[i + 1])
        metadata = 'Year %i-%i' % (i, i + 1)
        metric = metrics.HourglassMetric()
        displayDict = {'group': hourglassGroup, 'subgroup': 'Yearly', 'order': i}
        bundle = metricBundles.MetricBundle(metric, hourSlicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata)
        bundleList.append(bundle)

    # Histograms per filter for WFD only (generally used to produce merged histograms).
    summaryStats = standardStats
    prop = "WFD"
    for f in filters:
        # Set some per-proposal information.
        propCaption = ' for all WFD visits.'
        metadata = '%s band, WFD' % (f) + slicermetadata
        sqlconstraint = 'filter = "%s" and %s' % (f, wfdWhere)
        # Set up metrics and slicers for histograms.

        # Histogram the individual visit hour angle values.
        metric = metrics.CountMetric(col='HA', metricName='Hour Angle Histogram')
        histMerge = {'legendloc': 'upper right',
                'color': colors[f], 'label': '%s' % f,
                'xMin': -10., 'xMax': 10, "binsize": 0.1}
        plotDict = {"xMin": -10, "xMax": 10, "binsize": 0.1}
        displayDict = {'group': hourangleGroup, 'order': filtorder[f],
                       'caption': 'Histogram of the hour angle in %s band, %s' % (f, propCaption)}
        slicer = slicers.OneDSlicer(sliceColName='HA', binsize=0.1)
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata)
        mergedHistDict['hourAngle'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Histogram the individual visit distance to moon values.
        metric = metrics.CountMetric(col='moonDistance', metricName='Distance to Moon Histogram')
        histMerge = {'legendloc': 'upper right',
                     'color': colors[f], 'label': '%s' % f,
                     'xMin': 0., 'xMax': 180., "binsize": 5,
                     'xlabel': 'Distance to Moon (degrees)'}
        plotDict = {"xMin": 0, "xMax": 180, "binsize": 5}
        caption = 'Histogram of the distance between the field and the moon (in radians) '
        caption += 'in %s band, %s' % (f, propCaption)
        displayDict = {'group': distanceGroup, "subgroup": "Moon Distance",
                       'order': filtorder[f], 'caption': caption}
        slicer = slicers.OneDSlicer(sliceColName='moonDistance', binsize=5.)
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata)
        mergedHistDict['moonDistance'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

    # Slew histograms (time and distance).
    sqlconstraint = wfdWhere
    metadata = "WFD"
    metric = metrics.CountMetric(col='slewTime', metricName='Slew Time Histogram')
    plotDict = {'logScale': True, 'ylabel': 'Count',
                "xMin": -10, "xMax": 160, "binsize": 1}
    displayDict = {'group': slewGroup, 'subgroup': 'Slew Histograms',
                   'caption': 'Histogram of slew times for all visits.',
                   "order": 1}
    slicer = slicers.OneDSlicer(sliceColName='slewTime', binsize=1)
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName,
                                        metadata=metadata)
    bundleList.append(bundle)

    # now a zoomed in time count
    metric = metrics.CountMetric(col="slewTime", metricName="Slew Time Histogram (zoom)")
    plotDict = {"logScale": True, "ylabel": "Count",
                "xMin": 2, "xMax": 18, "binsize": 0.1}
    displayDict = {"group": slewGroup, "subgroup": "Slew Histograms",
                   "caption": "Zoom of histogram of slew times for all visits.",
                   "order": 2}
    slicer = slicers.OneDSlicer(sliceColName="slewTime", binsize=0.1)
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName,
                                        metadata=metadata)
    bundleList.append(bundle)

    # slew distance histogram
    metric = metrics.CountMetric(col='slewDistance', metricName='Slew Distance Histogram')
    plotDict = {'logScale': True, 'ylabel': 'Count',
                "xMin": -10, "xMax": 180, "binsize": 1}
    displayDict = {'group': slewGroup, 'subgroup': 'Slew Histograms',
                   'caption': 'Histogram of slew distances for all visits.',
                   "order": 3}
    slicer = slicers.OneDSlicer(sliceColName='slewDistance', binsize=3.)
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName,
                                        metadata=metadata)
    bundleList.append(bundle)

    # and a zoomed in slew distance histogram
    metric = metrics.CountMetric(col='slewDistance', metricName='Slew Distance Histogram (zoom)')
    plotDict = {'logScale': True, 'ylabel': 'Count',
                "xMin": 0, "xMax": 19, "binsize": 0.1}
    displayDict = {'group': slewGroup, 'subgroup': 'Slew Histograms',
                   'caption': 'Zoomed histogram of slew distances for all visits.',
                   "order": 4}
    slicer = slicers.OneDSlicer(sliceColName='slewDistance', binsize=0.1)
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName,
                                        metadata=metadata)
    bundleList.append(bundle)


    # Plots per night -- the number of visits and the open shutter time fraction.
    # nvisits per night
    slicer = slicers.OneDSlicer(sliceColName='night', binsize=1)
    metadata = 'Per night'
    sqlconstraint = ''
    summaryStats = allStats

    metric = metrics.CountMetric(col='observationStartMJD', metricName='NVisits')
    displayDict = {'group': summaryGroup, 'subgroup': '3: Obs Per Night',
                   'caption': 'Number of visits per night.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    # open shutter time per night
    metric = metrics.OpenShutterFractionMetric()
    caption = 'Open shutter fraction per night. This compares the on-sky image time against '
    caption += 'the on-sky time + slews/filter changes/readout, but does not include downtime '
    caption += 'due to weather.'
    displayDict = {'group': summaryGroup, 'subgroup': '3: Obs Per Night',
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    # filter stats per night
    metric = metrics.NChangesMetric(col='filter', metricName='Filter Changes')
    plotDict = {'ylabel': 'Number of Filter Changes'}
    displayDict = {'group': filterGroup, 'subgroup': 'Per Night',
                   'caption': 'Number of filter changes per night.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    metric = metrics.MinTimeBetweenStatesMetric(changeCol='filter')
    plotDict = {'yMin': 0, 'yMax': 120}
    displayDict = {'group': filterGroup, 'subgroup': 'Per Night',
                   'caption': 'Minimum time between filter changes, in minutes.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    cutoff = 10
    metric = metrics.NStateChangesFasterThanMetric(changeCol='filter', cutoff=cutoff)
    plotDict = {}
    caption = 'Number of filter changes, where the time between filter changes is shorter '
    caption += 'than %.1f minutes, per night.' % (cutoff)
    displayDict = {'group': filterGroup, 'subgroup': 'Per Night',
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    metric = metrics.NStateChangesFasterThanMetric(changeCol='filter', cutoff=20)
    plotDict = {}
    caption = 'Number of filter changes, where the time between filter changes '
    caption += 'is shorter than 20 minutes, per night.'
    displayDict = {'group': filterGroup, 'subgroup': 'Per Night',
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    metric = metrics.MaxStateChangesWithinMetric(changeCol='filter', timespan=10)
    plotDict = {}
    displayDict = {'group': filterGroup, 'subgroup': 'Per Night',
                   'caption': 'Max number of filter changes within a window of 10 minutes, per night.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    metric = metrics.MaxStateChangesWithinMetric(changeCol='filter', timespan=20)
    plotDict = {}
    displayDict = {'group': filterGroup, 'subgroup': 'Per Night',
                   'caption': 'Max number of filter changes within a window of 20 minutes, per night.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName, metadata=metadata,
                                        summaryMetrics=summaryStats)
    bundleList.append(bundle)

    # Unislicer (single number) metrics.
    slicer = slicers.UniSlicer()
    sqlconstraint = wfdWhere
    metadata = 'WFD only'
    order = 0

    metric = metrics.NChangesMetric(col='filter', metricName='Total Filter Changes')
    displayDict = {'group': filterGroup, 'subgroup': 'Whole Survey', 'order': order,
                   'caption': 'Total filter changes over survey'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    order += 1
    metric = metrics.MinTimeBetweenStatesMetric(changeCol='filter')
    displayDict = {'group': filterGroup, 'subgroup': 'Whole Survey', 'order': order,
                   'caption': 'Minimum time between filter changes, in minutes.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    order += 1
    metric = metrics.NStateChangesFasterThanMetric(changeCol='filter', cutoff=10)
    displayDict = {'group': filterGroup, 'subgroup': 'Whole Survey', 'order': order,
                   'caption': 'Number of filter changes faster than 10 minutes over the entire survey.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    order += 1
    metric = metrics.NStateChangesFasterThanMetric(changeCol='filter', cutoff=20)
    displayDict = {'group': filterGroup, 'subgroup': 'Whole Survey', 'order': order,
                   'caption': 'Number of filter changes faster than 20 minutes over the entire survey.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    order += 1
    metric = metrics.MaxStateChangesWithinMetric(changeCol='filter', timespan=10)
    caption = 'Max number of filter changes within a window of 10 minutes over the entire survey.'
    displayDict = {'group': filterGroup, 'subgroup': 'Whole Survey', 'order': order,
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    order += 1
    metric = metrics.MaxStateChangesWithinMetric(changeCol='filter', timespan=20)
    caption = 'Max number of filter changes within a window of 20 minutes over the entire survey.'
    displayDict = {'group': filterGroup, 'subgroup': 'Whole Survey', 'order': order,
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1

    # Calculate some basic summary info about run, per filter, per proposal and for all proposals.
    propOrder = 0
    prop = "WFD"
    slicer = slicers.UniSlicer()
    for f in filters + ['all']:
        if f != 'all':
            order = filtorder[f] * 100
            sqlconstraint = 'filter = "%s" and' % (f)
        else:
            order = 1000
            sqlconstraint = ''
        subgroup = 'WFD'
        sqlconstraint = sqlconstraint + ' %s' % (wfdWhere)
        metadata = '%s band, WFD' % (f)

        col = 'moonDistance'
        group = distanceGroup

        metric = metrics.MedianMetric(col=col)
        displayDict = {'group': group, 'subgroup': subgroup, 'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1
        metric = metrics.MeanMetric(col=col)
        displayDict = {'group': group, 'subgroup': subgroup, 'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1
        metric = metrics.RmsMetric(col=col)
        displayDict = {'group': group, 'subgroup': subgroup, 'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1
        metric = metrics.NoutliersNsigmaMetric(
            col=col, metricName='N(-3Sigma) %s' % (col), nSigma=-3.)
        displayDict = {'group': group, 'subgroup': subgroup, 'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1
        metric = metrics.NoutliersNsigmaMetric(col=col, metricName='N(+3Sigma) %s' % (col), nSigma=3.)
        displayDict = {'group': group, 'subgroup': subgroup, 'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1
        metric = metrics.CountMetric(col=col, metricName='Count %s' % (col))
        displayDict = {'group': group, 'subgroup': subgroup, 'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1
        metric = metrics.PercentileMetric(col=col, percentile=25)
        displayDict = {'group': group, 'subgroup': subgroup,
                       'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1
        metric = metrics.PercentileMetric(col=col, percentile=50)
        displayDict = {'group': group, 'subgroup': subgroup,
                       'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1
        metric = metrics.PercentileMetric(col=col, percentile=75)
        displayDict = {'group': group, 'subgroup': subgroup, 'order': order}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1

    # Calculate summary slew statistics.
    slicer = slicers.UniSlicer()
    sqlconstraint = ''
    metadata = 'All Visits'
    # Mean Slewtime
    metric = metrics.MeanMetric(col='slewTime')
    displayDict = {'group': slewGroup, 'subgroup': 'Summary', 'order': 1,
                   'caption': 'Mean slew time in seconds.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)
    # Median Slewtime
    metric = metrics.MedianMetric(col='slewTime')
    displayDict = {'group': slewGroup, 'subgroup': 'Summary', 'order': 2,
                   'caption': 'Median slew time in seconds.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)
    # Mean exposure time
    metric = metrics.MeanMetric(col='visitExposureTime')
    displayDict = {'group': slewGroup, 'subgroup': 'Summary', 'order': 3,
                   'caption': 'Mean visit on-sky time, in seconds.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)
    # Mean visit time
    metric = metrics.MeanMetric(col='visitTime')
    displayDict = {'group': slewGroup, 'subgroup': 'Summary', 'order': 4,
                   'caption':
                   'Mean total visit time (including readout and shutter), in seconds.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    # Count total nvisits
    sqlconstraint = ''
    slicer = slicers.UniSlicer()
    metadata = 'All Visits'

    metric = metrics.CountMetric(col='observationStartMJD', metricName='NVisits')
    summaryMetrics = [metrics.IdentityMetric(metricName='Count')]
    displayDict = {'group': summaryGroup, 'subgroup': '0: NVisits', 'order': 0}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        summaryMetrics=summaryMetrics,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    # and total exposure time
    sqlconstraint = ''
    slicer = slicers.UniSlicer()
    metadata = 'All Visits'

    metric = NormSumMetric(col='visitExposureTime', normValue=3600*24,
                           metricName='Exp. Time (days)')
    summaryMetrics = [metrics.IdentityMetric(metricName='days')]
    displayDict = {'group': summaryGroup, 'subgroup': '1: Exp. Time', 'order': 1}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        summaryMetrics=summaryMetrics,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)



    # Count the total exposure time per proposal, for all proposals
    order = 2
    slicer = slicers.UniSlicer()
    for propid in propids:
        sqlconstraint = '%s = %s' % (propCol, propid)
        metadata = '%s' % (propids[propid])

        metric = metrics.CountMetric(col='observationStartMJD', metricName='NVisits')
        summaryMetrics = [metrics.IdentityMetric(metricName='Count')]
        displayDict = {'group': summaryGroup, 'subgroup': '0: NVisits', 'order': order,
                       'caption': 'Number of visits for %s proposal.'
                       % (propids[propid])}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            summaryMetrics=summaryMetrics,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)

        order += 1

        metric = NormSumMetric(col="visitExposureTime", normValue=3600*24,
                               metricName="Exp. Time (days)")
        summaryMetrics = [metrics.IdentityMetric(metricName='days')]
        displayDict = {'group': summaryGroup, 'subgroup': '1: Exp. Time', 'order': order,
                       'caption': 'Total exposure time for %s proposal (days).'
                       % (propids[propid])}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                            summaryMetrics=summaryMetrics,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata)
        bundleList.append(bundle)


        order += 1

    # Count total number of nights
    metadata = "All Visits"
    sqlconstraint = ''
    metric = metrics.CountUniqueMetric(col='night', metricName='Nights with observations')
    summaryMetrics = [metrics.IdentityMetric(metricName='(days)')]
    displayDict = {'group': summaryGroup, 'subgroup': '2: On-sky Time', 'order': 1}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        summaryMetrics=summaryMetrics,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    metric = metrics.FullRangeMetric(col='night', metricName='Total nights in survey')
    summaryMetrics = [metrics.ZeropointMetric(zp=1, metricName='(days)')]
    displayDict = {'group': summaryGroup, 'subgroup': '2: On-sky Time', 'order': 0}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        summaryMetrics=summaryMetrics,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    metric = metrics.TeffMetric(metricName='Total effective time of survey')
    summaryMetrics = [metrics.NormalizeMetric(normVal=24.0 * 60.0 * 60.0, metricName='(days)')]
    displayDict = {'group': summaryGroup, 'subgroup': '2: On-sky Time', 'order': 3}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        summaryMetrics=summaryMetrics,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    metric = metrics.TeffMetric(metricName='Normalized total effective time of survey', normed=True)
    summaryMetrics = [metrics.IdentityMetric(metricName='(fraction)')]
    displayDict = {'group': summaryGroup, 'subgroup': '2: On-sky Time', 'order': 2}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        summaryMetrics=summaryMetrics,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)

    # Check the Alt-Az pointing history
    slicer = slicers.HealpixSlicer(nside=nside, latCol='zenithDistance', lonCol='azimuth', useCache=False)
    metric = metrics.CountMetric('observationStartMJD', metricName='NVisits Alt/Az')
    plotDict = {'rot': (0, 90, 0)}
    plotFunc = plots.HealpixSkyMap()
    for f in filters:
        sqlconstraint = 'filter = "%s"' % (f)
        displayDict = {'group': hourangleGroup, 'order': filtorder[f],
                       'caption':
                       'Pointing History on the alt-az sky (zenith center) for filter %s' % f}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            runName=runName,
                                            plotFuncs=[plotFunc], displayDict=displayDict)
        bundleList.append(bundle)
    displayDict = {'group': hourangleGroup, 'subgroup': 'All Filters',
                   'caption':
                   'Pointing History on the alt-az sky (zenith center), all filters'}
    bundle = metricBundles.MetricBundle(metric, slicer, '', plotDict=plotDict, runName=runName,
                                        plotFuncs=[plotFunc], displayDict=displayDict)
    bundleList.append(bundle)

    # Solar elongation
    sqls = ['filter = "%s"' % f for f in filters]
    orders = [filtorder[f] for f in filters]
    sqls.append('')
    orders.append(0)
    for sql, order in zip(sqls, orders):
        plotFuncs = [plots.HealpixSkyMap(), plots.HealpixHistogram()]
        displayDict = {'group': distanceGroup, 'subgroup': 'Solar Elongation',
                       'caption': 'Median solar elongation in degrees', 'order': order}
        metric = metrics.MedianMetric('solarElong')
        slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
        bundle = metricBundles.MetricBundle(metric, slicer, sql, displayDict=displayDict,
                                            plotFuncs=plotFuncs, runName=runName)
        bundleList.append(bundle)


    return (metricBundles.makeBundlesDictFromList(bundleList), mergedHistDict)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Python script to run MAF with the scheduler validation metrics')
    parser.add_argument('dbFile', type=str, default=None, help="full file path to the opsim sqlite file")

    parser.add_argument("--outDir", type=str, default='./Out', help='Output directory for MAF outputs.')

    parser.add_argument('--benchmark', type=str, default='design',
                        help="Can be 'design' or 'requested'")

    parser.add_argument('--plotOnly', dest='plotOnly', action='store_true',
                        default=False, help="Reload the metric values and re-plot them.")

    parser.add_argument('--skipSlew', dest='skipSlew', action='store_true',
                        default=False, help='Skip calculation of slew statistics')
    parser.add_argument("--runName", type=str, help="Name of the run", default=None)

    parser.set_defaults()
    args, extras = parser.parse_known_args()

    resultsDb = db.ResultsDb(outDir=args.outDir)
    opsdb = db.OpsimDatabaseV4(args.dbFile)

    (bundleDict, mergedHistDict) = makeBundleList(args.dbFile, runName=args.runName, benchmark=args.benchmark)

    group = metricBundles.MetricBundleGroup(bundleDict, opsdb, outDir=args.outDir, resultsDb=resultsDb)
    if args.plotOnly:
        group.readAll()
    else:
        group.runAll()
    group.plotAll()

    for key in mergedHistDict:
        if len(mergedHistDict[key].bundleList) > 0:
            mergedHistDict[key].percentileLegend()
            mergedHistDict[key].incrementPlotOrder()
            mergedHistDict[key].plot(outDir=args.outDir, resultsDb=resultsDb, closeFigs=True)
        else:
            warnings.warn('Empty bundleList for %s, skipping merged histogram' % key)

    utils.writeConfigs(opsdb, args.outDir)
    opsdb.close()
