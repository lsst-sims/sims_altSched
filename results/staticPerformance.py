#!/usr/bin/env python
from __future__ import print_function
from builtins import zip
import os
import argparse
import warnings
# Set matplotlib backend (to create plots where DISPLAY is not set).
import matplotlib
matplotlib.use('Agg')
import numpy as np
import healpy as hp
import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import copy

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.utils as utils


def makeBundleList(dbFile, runName=None, nside=64, benchmark='design',
                   lonCol='fieldRA', latCol='fieldDec'):
    """
    make a list of metricBundle objects to look at the scientific performance
    of an opsim run.
    """

    # List to hold everything we're going to make
    bundleList = []

    # Connect to the databse
    opsimdb = db.OpsimDatabaseV4(dbFile)
    if runName is None:
        runName = os.path.basename(dbFile).replace('_sqlite.db', '')

    # Fetch the proposal ID values from the database
    propids, propTags = opsimdb.fetchPropInfo()

    # Fetch the telescope location from config
    lat, lon, height = opsimdb.fetchLatLonHeight()

    # Construct a WFD SQL where clause so multiple propIDs can query by WFD:
    wfdWhere = opsimdb.createSQLWhere('WFD', propTags)
    print('#FYI: WFD "where" clause: %s' % (wfdWhere))
    ddWhere = opsimdb.createSQLWhere('DD', propTags)
    print('#FYI: DD "where" clause: %s' % (ddWhere))

    # Set up benchmark values, scaled to length of opsim run.
    runLength = opsimdb.fetchRunLength()
    benchmarkSeeing = 'FWHMeff'
    if benchmark == 'requested':
        # Fetch design values for seeing/skybrightness/single visit depth.
        benchmarkVals = utils.scaleBenchmarks(runLength, benchmark='design')
        # Update nvisits with requested visits from config files.
        benchmarkVals['nvisits'] = opsimdb.fetchRequestedNvisits(propId=propTags['WFD'])
        # Calculate expected coadded depth.
        benchmarkVals['coaddedDepth'] = utils.calcCoaddedDepth(benchmarkVals['nvisits'],
                                                               benchmarkVals['singleVisitDepth'])
    elif (benchmark == 'stretch') or (benchmark == 'design'):
        # Calculate benchmarks for stretch or design.
        benchmarkVals = utils.scaleBenchmarks(runLength, benchmark=benchmark)
        benchmarkVals['coaddedDepth'] = utils.calcCoaddedDepth(benchmarkVals['nvisits'],
                                                               benchmarkVals['singleVisitDepth'])
    else:
        raise ValueError('Could not recognize benchmark value %s, use design, stretch or requested.'
                         % (benchmark))
    # Check that nvisits is not set to zero (for very short run length).
    for f in benchmarkVals['nvisits']:
        if benchmarkVals['nvisits'][f] == 0:
            print('Updating benchmark nvisits value in %s to be nonzero' % (f))
            benchmarkVals['nvisits'][f] = 1

    # Set values for min/max range of nvisits for All/WFD and DD plots. These are somewhat arbitrary.
    nvisitsRange = {}
    nvisitsRange['all'] = {'u': [30, 110], 'g': [50, 160], 'r': [100, 310],
                           'i': [90, 270], 'z': [100, 330], 'y': [50, 300]}

    # Scale these ranges for the runLength.
    scale = runLength / 10.0
    for prop in nvisitsRange:
        for f in nvisitsRange[prop]:
            for i in [0, 1]:
                nvisitsRange[prop][f][i] = int(np.floor(nvisitsRange[prop][f][i] * scale))

    # Filter list, and map of colors (for plots) to filters.
    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    colors = {'u': 'cyan', 'g': 'g', 'r': 'y', 'i': 'r', 'z': 'm', 'y': 'k'}
    filtorder = {'u': 1, 'g': 2, 'r': 3, 'i': 4, 'z': 5, 'y': 6}

    # Easy way to run through all fi

    # Set up a list of common summary stats
    commonSummary = [metrics.MeanMetric(), metrics.RobustRmsMetric(), metrics.MedianMetric(),
                     metrics.PercentileMetric(metricName='25th%ile', percentile=25),
                     metrics.PercentileMetric(metricName='75th%ile', percentile=75),
                     metrics.MinMetric(), metrics.MaxMetric()]

    standardStats = [metrics.MeanMetric(),
                     metrics.RmsMetric(), metrics.MedianMetric(), metrics.CountMetric(),
                     metrics.NoutliersNsigmaMetric(metricName='N(+3Sigma)', nSigma=3),
                     metrics.NoutliersNsigmaMetric(metricName='N(-3Sigma)', nSigma=-3.)]

    allStats = copy.deepcopy(standardStats)
    allStats.extend(commonSummary)

    # Set up some 'group' labels
    summaryGroup = "A: Summary"
    fOGroup = 'B: fO'
    depthGroup = 'C: Depth per filter'
    airmassGroup = 'D: Airmass distribution'
    seeingGroup = 'E: Seeing distribution'
    skyBrightnessGroup = "F: Sky brightness distribution"
    singleDepthGroup = "G: Single Visit Depths"

    # Set up an object to track the metricBundles that we want to combine into merged plots.
    mergedHistDict = {}

    # Set the histogram merge function.
    mergeFunc = plots.HealpixHistogram()

    healPlots = ['coaddM5', 'normTEff', 'medianSingleDepth',
                 'medianSkyBrightness', 'medianSeeing', 'seeingAboveLimit',
                 'medianAirmass', 'medianNormAirmass'] 
    oneDPlots = ['singleDepth', 'skyBrightness', 'seeing', 
                 'airmass', 'normAirmass']
    
    # set the mergedHistDict valus to healpix and oneD plot bundles depending on
    # which plot is being made (median<...> measure median values to a sky pixel;
    # <...> without the median are just a histogram over all visits)
    for plotName in healPlots:
        mergedHistDict[plotName] = plots.PlotBundle(plotFunc=plots.HealpixHistogram())
    for plotName in oneDPlots:
        mergedHistDict[plotName] = plots.PlotBundle(plotFunc=plots.OneDBinnedData())

    ##
    # Calculate the fO metrics for all proposals and WFD only.
    order = 0
    prop = "WFD only"
    metadata = 'WFD only'
    sqlconstraint = '%s' % (wfdWhere)
    # Configure the count metric which is what is used for f0 slicer.
    m1 = metrics.CountMetric(col='observationStartMJD', metricName='fO')
    plotDict = {'xlabel': 'Number of Visits', 'Asky': benchmarkVals['Area'],
                'Nvisit': benchmarkVals['nvisitsTotal'], 'xMin': 0, 'xMax': 1500}
    summaryMetrics = [metrics.fOArea(nside=nside, norm=False, metricName='fOArea: Nvisits (#)',
                                     Asky=benchmarkVals['Area'], Nvisit=benchmarkVals['nvisitsTotal']),
                      metrics.fOArea(nside=nside, norm=True, metricName='fOArea: Nvisits/benchmark',
                                     Asky=benchmarkVals['Area'], Nvisit=benchmarkVals['nvisitsTotal']),
                      metrics.fONv(nside=nside, norm=False, metricName='fONv: Area (sqdeg)',
                                   Asky=benchmarkVals['Area'], Nvisit=benchmarkVals['nvisitsTotal']),
                      metrics.fONv(nside=nside, norm=True, metricName='fONv: Area/benchmark',
                                   Asky=benchmarkVals['Area'], Nvisit=benchmarkVals['nvisitsTotal'])]
    caption = 'The FO metric evaluates the overall efficiency of observing. '
    caption += ('fOArea: Nvisits = %.1f sq degrees receive at least this many visits out of %d. '
                % (benchmarkVals['Area'], benchmarkVals['nvisitsTotal']))
    caption += ('fONv: Area = this many square degrees out of %.1f receive at least %d visits.'
                % (benchmarkVals['Area'], benchmarkVals['nvisitsTotal']))
    displayDict = {'group': fOGroup, 'subgroup': 'F0', 'displayOrder': order, 'caption': caption}
    order += 1
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)

    bundle = metricBundles.MetricBundle(m1, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryMetrics,
                                        plotFuncs=[plots.FOPlot()],
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)

    order += 1

    ##
    # Depth metrics.
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
    for f in filters:
        propCaption = '%s band, all proposals' % (f)
        sqlconstraint = 'filter = "%s"' % (f)
        metadata = '%s band' % (f) 
        # Coadded depth.
        metric = metrics.Coaddm5Metric()
        plotDict = {'zp': benchmarkVals['coaddedDepth'][f],
                    'xlabel': 'coadded m5 - %.1f' % benchmarkVals['coaddedDepth'][f],
                    'xMin':     -1.5, 'xMax':     0.5, "binsize": 0.02,
                    "colorMin": -1.5, "colorMax": 0.5}
        summaryStats = allStats
        histMerge = {'legendloc': 'upper right', 'color': colors[f], 'label': '%s' % f, 'binsize': .02,
                     'xlabel': 'coadded m5 - benchmark value',
                     "xMin": -1.5, "xMax": 0.5, "binsize": 0.02}
        caption = ('Coadded depth in filter %s, with %s value subtracted (%.1f), %s. '
                   % (f, benchmark, benchmarkVals['coaddedDepth'][f], propCaption))
        caption += 'More positive numbers indicate fainter limiting magnitudes.'
        displayDict = {'group': depthGroup, 'subgroup': 'Coadded Depth',
                       'order': filtorder[f], 'caption': caption}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata,
                                            summaryMetrics=summaryStats)
        mergedHistDict['coaddM5'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)
        # Effective time.
        metric = metrics.TeffMetric(metricName='Normalized Effective Time', normed=True,
                                    fiducialDepth=benchmarkVals['singleVisitDepth'])
        plotDict = {'xMin':     0, 'xMax':     1.1, "binsize": 0.01,
                    "colorMin": 0, "colorMax": 1.1}
        summaryStats = allStats
        histMerge = {'legendLoc': 'upper right', 'color': colors[f], 'label': '%s' % f, 'binsize': 0.02}
        caption = ('"Time Effective" in filter %s, calculated with fiducial single-visit depth of %s mag. '
                   % (f, benchmarkVals['singleVisitDepth'][f]))
        caption += 'Normalized by the fiducial time effective, if every observation was at '
        caption += 'the fiducial depth.'
        displayDict = {'group': depthGroup, 'subgroup': 'Time Eff.',
                       'order': filtorder[f], 'caption': caption}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata,
                                            summaryMetrics=summaryStats)
        mergedHistDict['normTEff'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

    # Good seeing in r/i band metrics, including in first/second years.
    order = 0
    for tcolor, tlabel, timespan in zip(['k', 'g', 'r'], ['10 years', '1 year', '2 years'],
                                        ['', ' and night<=365', ' and night<=730']):
        order += 1
        for f in (['r', 'i']):
            sqlconstraint = 'filter = "%s" %s' % (f, timespan)
            propCaption = '%s band, all proposals, over %s.' % (f, tlabel)
            metadata = '%s band, %s' % (f, tlabel)

            seeing_limit = 0.7

            metric = metrics.FracAboveMetric(col="seeingFwhmEff", cutoff=seeing_limit)
            summaryStats = allStats
            plotDict = {'xMin': 0, 'xMax': 1.1, "binsize": 0.01,
                        'colorMin': 0, 'colorMax': 1.1, 'color': tcolor}
            displayDict = {'group': seeingGroup, 'subgroup': 'Good seeing fraction',
                           'order': filtorder[f] * 100 + order,
                           'caption': 'Fraction of total images with FWHMEff worse than %.1f, in %s'
                           % (seeing_limit, propCaption)}
            histMerge = {'color': tcolor, 'label': '%s %s' % (f, tlabel),
                         'binsize': 0.05, 'legendloc': 'upper right'}
            bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                                displayDict=displayDict, runName=runName, metadata=metadata,
                                                summaryMetrics=summaryStats)
            mergedHistDict['seeingAboveLimit'].addBundle(bundle, plotDict=histMerge)
            bundleList.append(bundle)

    prop = "WFD"
    for f in filters:
        propCaption = 'for all WFD proposals'
        metadata = '%s band, WFD' % (f)
        sqlconstraint = 'filter = "%s" and %s' % (f, wfdWhere)
        nvisitsMin = nvisitsRange['all'][f][0]
        nvisitsMax = nvisitsRange['all'][f][1]
        mag_zp = benchmarkVals['coaddedDepth'][f]

        slicer = slicers.HealpixSlicer()

        # Calculate the median individual visit five sigma limiting magnitude
        # (individual image depth).
        metric = metrics.MedianMetric(col='fiveSigmaDepth')
        summaryStats = standardStats
        plotDict = {"xMin":     20, "xMax":     26,
                    "colorMin": 20, "colorMax": 26,
                    "binsize": 0.05}
        displayDict = {'group': singleDepthGroup, 'order': filtorder[f],
                       'caption': 'Median single visit depth in filter %s, %s.' % (f, propCaption)}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata,
                                            summaryMetrics=summaryStats)
        histMerge = {'color': colors[f], 'label': '%s' % (f),
                     'xlabel': 'Median 5-sigma depth (mags)',
                     'binsize': .05, 'legendloc': 'upper right',
                     'xMin': 20, 'xMax': 26}
        mergedHistDict['medianSingleDepth'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Calculate the median individual visit sky brightness (normalized to a benchmark).
        metric = metrics.MedianMetric(col='skyBrightness')
        plotDict = {'zp': benchmarkVals['skybrightness'][f],
                    'xlabel': 'Skybrightness - %.2f' % (benchmarkVals['skybrightness'][f]),
                    'xMin':     -1.5, 'xMax':     1.5, "binsize": 0.02,
                    'colorMin': -1.5, 'colorMax': 1.5,
                    'cmap': cm.RdBu}
        caption = 'Median Sky Brightness in filter %s ' % f
        caption += 'with expected zeropoint (%.2f) subtracted, ' % (benchmarkVals['skybrightness'][f])
        caption += '%s. ' % (propCaption)
        caption += 'Fainter sky brightness values are more positive numbers.'
        displayDict = {'group': skyBrightnessGroup, 'order': filtorder[f],
                       'caption': caption}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata,
                                            summaryMetrics=summaryStats)
        histMerge = {'zp': benchmarkVals['skybrightness'][f], 'color': colors[f],
                     'label': '%s' % (f),
                     'xMin': -1.5, 'xMax': 1.5,
                     'binsize': .02, 
                     'xlabel': 'Skybrightness - benchmark',
                     'legendloc': 'upper right'}
        mergedHistDict['medianSkyBrightness'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Calculate the median delivered seeing.
        metric = metrics.MedianMetric(col="seeingFwhmEff")
        plotDict = {'normVal': benchmarkVals[benchmarkSeeing][f],
                    'xlabel': 'Median FWHMEff/(Expected FWHMEff %.2f)'
                    % (benchmarkVals[benchmarkSeeing][f]),
                    'xMin':     0.8, 'xMax':     1.6,
                    'colorMin': 0.8, 'colorMax': 1.6,
                    'binsize': 0.02, 'cmap': cm.RdBu_r}
        caption = 'Median Seeing in filter %s ' % f
        caption += 'divided by expected value (%.2f), %s. ' % (benchmarkVals[benchmarkSeeing][f],
                                                               propCaption)
        caption += 'Seeing is seeingFwhmEff column.'
        displayDict = {'group': seeingGroup, 'order': filtorder[f], 'caption': caption}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata,
                                            summaryMetrics=summaryStats)
        histMerge = {'color': colors[f], 'label': '%s' % (f),
                     'xlabel': 'Seeing/benchmark seeing',
                     'xMin': 0.8, 'xMax': 1.6,
                     'binsize': .02,
                     'legendloc': 'upper right'}
        mergedHistDict['medianSeeing'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Calculate the median airmass.
        metric = metrics.MedianMetric(col='airmass')
        plotDict = {'xMin':     1, 'xMax':     2,
                    'colorMin': 1, 'colorMax': 2,
                    'binsize': 0.01}
        displayDict = {'group': airmassGroup, 'order': filtorder[f],
                       'caption': 'Median airmass in filter %s, %s.' % (f, propCaption)}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata,
                                            summaryMetrics=summaryStats)
        histMerge = {'color': colors[f], 'label': '%s' % (f),
                     'xlabel': 'Median Airmass', 'binsize': .01,
                     'xMin': 1, 'xMax': 2,
                     'legendloc': 'upper right'}
        mergedHistDict['medianAirmass'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Calculate the median normalized airmass.
        metric = metrics.MedianMetric(col='normairmass')
        plotDict = {'xMin':     1, 'xMax':     1.5,
                    'colorMin': 1, 'colorMax': 1.5,
                    'binsize': 0.01}
        caption = 'Median normalized airmass (airmass divided by the minimum airmass '
        caption += 'a field could reach) in filter %s, %s.' % (f, propCaption)
        displayDict = {'group': airmassGroup, 'order': filtorder[f],
                       'caption': caption}
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName,
                                            metadata=metadata,
                                            summaryMetrics=summaryStats)
        histMerge = {'color': colors[f], 'label': '%s' % (f),
                     'xlabel': 'Median Normalized Airmass',
                     'xMin': 1, 'xMax': 1.5, 'binsize': .01,
                     'legendloc': 'upper right'}
        mergedHistDict['medianNormAirmass'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

    # Histograms per filter for WFD only (generally used to produce merged histograms).
    summaryStats = standardStats
    prop = "WFD"
    for f in filters:
        # Set some per-proposal information.
        propCaption = ' for all WFD visits.'
        metadata = '%s band, WFD' % (f)
        sqlconstraint = 'filter = "%s" and %s' % (f, wfdWhere)

        # Histogram the individual visit five sigma limiting magnitude (individual image depth).
        metric = metrics.CountMetric(col='fiveSigmaDepth', metricName='Single Visit Depth Histogram')
        histMerge = {'legendloc': 'upper right', 'color': colors[f], 'label': '%s' % f,
                     "xMin": 20, "xMax": 26, "binsize": 0.05}
        caption = 'Histogram of the single visit depth in %s band, %s.' % (f, propCaption)
        plotDict = {"xMin": 20, "xMax": 26, "binsize": 0.05}
        displayDict = {'group': singleDepthGroup, 'order': filtorder[f],
                       'caption': caption}
        slicer = slicers.OneDSlicer(sliceColName='fiveSigmaDepth', binsize=0.05)
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata)
        mergedHistDict['singleDepth'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Histogram the individual visit sky brightness.
        metric = metrics.CountMetric(col='skyBrightness', metricName='Sky Brightness Histogram')
        histMerge = {'legendloc': 'upper right',
                     'color': colors[f], 'label': '%s' % f,
                     "xMin": 16, "xMax": 24, "binsize": 0.05}
        plotDict = {"xMin": 16, "xMax": 24, "binsize": 0.05}
        displayDict = {'group': skyBrightnessGroup, 'order': filtorder[f],
                       'caption': 'Histogram of the sky brightness in %s band, %s.' % (f, propCaption)}
        slicer = slicers.OneDSlicer(sliceColName='skyBrightness', binsize=0.05,
                                    binMin=16, binMax=24)
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata)
        mergedHistDict['skyBrightness'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Histogram the individual visit seeing.
        metric = metrics.CountMetric(col="seeingFwhmEff", metricName='Seeing Histogram')
        histMerge = {'legendloc': 'upper right',
                     'color': colors[f], 'label': '%s' % f,
                     "xMin": 0.4, "xMax": 2.5, "binsize": 0.02}
        caption = 'Histogram of the seeing in %s band, %s. ' % (f, propCaption)
        caption += 'Seeing is fwhmEff column.'
        displayDict = {'group': seeingGroup, 'order': filtorder[f],
                       'caption': caption}
        plotDict = {"xMin": 0.4, "xMax": 2.5, "binsize": 0.02}
        slicer = slicers.OneDSlicer(sliceColName="seeingFwhmEff", binsize=0.02)
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata)
        mergedHistDict['seeing'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Histogram the individual visit airmass values.
        metric = metrics.CountMetric(col='airmass', metricName='Airmass Histogram')
        histMerge = {'legendloc': 'upper right',
                     'color': colors[f], 'label': '%s' % f, 'xMin': 1.0, 'xMax': 2.0}
        plotDict = {"xMin": 1, "xMax": 2, "binsize": 0.01}
        displayDict = {'group': airmassGroup, 'order': filtorder[f],
                       'caption': 'Histogram of the airmass in %s band, %s' % (f, propCaption)}
        slicer = slicers.OneDSlicer(sliceColName='airmass', binsize=0.01)
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata)
        mergedHistDict['airmass'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

        # Histogram the individual visit normalized airmass values.
        metric = metrics.CountMetric(col='normairmass', metricName='Normalized Airmass Histogram')
        histMerge = {'legendloc': 'upper right',
                     'color': colors[f], 'label': '%s' % f,
                     'xMin': 1.0, 'xMax': 1.4, "binsize": 0.01}
        plotDict = {"xMin": 1, "xMax": 1.4}
        displayDict = {'group': airmassGroup, 'order': filtorder[f],
                       'caption': 'Histogram of the normalized airmass in %s band, %s' % (f, propCaption)}
        slicer = slicers.OneDSlicer(sliceColName='normairmass', binsize=0.01)
        bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName, metadata=metadata)
        mergedHistDict['normAirmass'].addBundle(bundle, plotDict=histMerge)
        bundleList.append(bundle)

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

        cols = ["seeingFwhmEff", 'skyBrightness', 'airmass',
                'fiveSigmaDepth', 'normairmass']
        groups = [seeingGroup, skyBrightnessGroup, airmassGroup,
                  singleDepthGroup, airmassGroup]
        for col, group in zip(cols, groups):
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

    # total exposure time of survey
    sqlconstraint = ''
    slicer = slicers.UniSlicer()
    metadata = 'All Visits'

    metric = metrics.SumMetric(col='visitExposureTime', metricName='Total Exposure Time (days)')
    summaryMetrics = [metrics.NormalizeMetric(normVal=24*3600, metricName='(days)')]
    displayDict = {'group': summaryGroup, 'subgroup': '1: NVisits', 'order': 0}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        summaryMetrics=summaryMetrics,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)


    # total T_eff of the survey
    metric = metrics.TeffMetric(metricName='Total effective time of survey (days)')
    summaryMetrics = [metrics.NormalizeMetric(normVal=24.0 * 60.0 * 60.0, metricName='(days)')]
    displayDict = {'group': summaryGroup, 'subgroup': '2: On-sky Time', 'order': 3}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint,
                                        summaryMetrics=summaryMetrics,
                                        displayDict=displayDict, runName=runName, metadata=metadata)
    bundleList.append(bundle)


    return (metricBundles.makeBundlesDictFromList(bundleList), mergedHistDict)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Python script to run the science performance metrics.')
    parser.add_argument('dbFile', type=str, default=None, help="full file path to the opsim sqlite file")
    parser.add_argument("--outDir", type=str, default='./Out', help='Output directory for MAF outputs.' +
                        ' Default "Out"')
    parser.add_argument("--nside", type=int, default=64,
                        help="Resolution to run Healpix grid at (must be 2^x). Default 64.")
    parser.add_argument("--lonCol", type=str, default='fieldRA',
                        help="Column to use for RA values (can be a stacker dither column)." +
                        " Default=fieldRA.")
    parser.add_argument("--latCol", type=str, default='fieldDec',
                        help="Column to use for Dec values (can be a stacker dither column)." +
                        " Default=fieldDec.")
    parser.add_argument('--benchmark', type=str, default='design',
                        help="Can be 'design' or 'requested'")
    parser.add_argument('--plotOnly', dest='plotOnly', action='store_true',
                        default=False, help="Reload the metric values from disk and re-plot them.")
    parser.add_argument('--skipNoSave', dest='runNoSave', action='store_false', default=True,
                        help="Skip the metrics that do not get saved as npz files.")
    parser.add_argument("--runName", type=str, help="Name of the run")
    parser.set_defaults()
    args, extras = parser.parse_known_args()

    # Build metric bundles.

    (bundleDict, mergedHistDict) = makeBundleList(args.dbFile, nside=args.nside,
                                                  lonCol=args.lonCol, latCol=args.latCol,
                                                  benchmark=args.benchmark, runName=args.runName)

    # Set up / connect to resultsDb.
    resultsDb = db.ResultsDb(outDir=args.outDir)
    # Connect to opsimdb.
    opsdb = db.OpsimDatabaseV4(args.dbFile)

    # Set up metricBundleGroup.
    group = metricBundles.MetricBundleGroup(bundleDict, opsdb,
                                            outDir=args.outDir, resultsDb=resultsDb)
    # Read or run to get metric values.
    if args.plotOnly:
        group.readAll()
    else:
        group.runAll()
    # Make plots.
    group.plotAll()
    # Make merged plots.
    for key in mergedHistDict:
        if len(mergedHistDict[key].bundleList) > 0:
            mergedHistDict[key].incrementPlotOrder()
            mergedHistDict[key].plot(outDir=args.outDir, resultsDb=resultsDb, closeFigs=True)
        else:
            warnings.warn('Empty bundleList for %s, skipping merged histogram' % key)
    # Get config info and write to disk.
    utils.writeConfigs(opsdb, args.outDir)
    opsdb.close()

    print("Finished sciencePerformance metric calculations.")
