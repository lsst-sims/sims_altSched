from __future__ import print_function, division

import matplotlib
matplotlib.use("agg")
import argparse
import numpy as np

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.plots as plots
import lsst.sims.maf.utils as utils

from timeMetrics import LonelinessMetric, IntranightColorMetric, SNLotsMetric

def makeBundleList(dbFile, runName, nside=64):
    bundleList = []

    opsimdb = db.OpsimDatabaseV4(dbFile)

    # Fetch the proposal ID values from the database
    propids, propTags = opsimdb.fetchPropInfo()

    # Fetch the telescope location from config
    lat, lon, height = opsimdb.fetchLatLonHeight()

    # Construct a WFD SQL where clause so multiple propIDs can query by WFD:
    wfdWhere = opsimdb.createSQLWhere('WFD', propTags)
    print('#FYI: WFD "where" clause: %s' % (wfdWhere))
    ddWhere = opsimdb.createSQLWhere('DD', propTags)
    print('#FYI: DD "where" clause: %s' % (ddWhere))

    # Filter list, and map of colors (for plots) to filters.
    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    colors = {'u': 'cyan', 'g': 'g', 'r': 'y', 'i': 'r', 'z': 'm', 'y': 'k'}
    filtorder = {'u': 1, 'g': 2, 'r': 3, 'i': 4, 'z': 5, 'y': 6}
 
    # Set up a list of common summary stats
    commonSummary = [metrics.MeanMetric(), metrics.RobustRmsMetric(), metrics.MedianMetric(),
                     metrics.PercentileMetric(metricName='25th%ile', percentile=25),
                     metrics.PercentileMetric(metricName='75th%ile', percentile=75),
                     metrics.MinMetric(), metrics.MaxMetric()]
    allStats = commonSummary

    parallaxGroup = "A: Parallax Metrics"
    properMotionGroup = "B: Proper Motion Metrics"
    intraNightGroup = "C: Intra-Night Visit Statistics"
    interNightGroup = "D: Inter-Night Visit Statistics"
    snGroup = "E: Type Ia Supernovae"
    neoGroup = "F: NEO Metrics"

    latCol = "fieldDec"
    lonCol = "fieldRA"
    geomSeeingCol = "seeingFwhmGeom"

    ##
    # Trigonometric parallax and proper motion @ r=20 and r=24
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
    sqlconstraint = ''
    metadata = "All Visits"
    order = 0
    metric = metrics.ParallaxMetric(metricName='Parallax 20', rmag=20,
                                    seeingCol=geomSeeingCol)
    summaryStats = allStats
    plotDict = {'cbarFormat': '%.1f',
                'xMin':     0, 'xMax':     2,
                'colorMin': 0, 'colorMax': 2,
                'binsize': 0.05}
    displayDict = {'group': parallaxGroup, 'subgroup': 'Parallax', 'order': order,
                   'caption': 'Parallax precision at r=20. (without refraction).'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1

    metric = metrics.ParallaxMetric(metricName='Parallax 24', rmag=24,
                                    seeingCol=geomSeeingCol)
    plotDict = {'cbarFormat': '%.1f',
                'xMin':     2, 'xMax':     8,
                'colorMin': 2, 'colorMax': 8,
                'binsize': 0.05}
    displayDict = {'group': parallaxGroup, 'subgroup': 'Parallax', 'order': order,
                   'caption': 'Parallax precision at r=24. (without refraction).'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1
    metric = metrics.ParallaxMetric(metricName='Parallax Normed', rmag=24, normalize=True,
                                    seeingCol=geomSeeingCol)
    plotDict = {'xMin':     0, 'xMax':     1, "binsize": 0.01,
                "colorMin": 0, "colorMax": 1}
    displayDict = {'group': parallaxGroup, 'subgroup': 'Parallax', 'order': order,
                   'caption':
                   'Normalized parallax (normalized to optimum observation cadence, 1=optimal).'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1
    metric = metrics.ParallaxCoverageMetric(metricName='Parallax Coverage 20',
                                            rmag=20, seeingCol=geomSeeingCol)
    plotDict = {"xMin":     0, "xMax":     1, "binsize": 0.01,
                "colorMin": 0, "colorMax": 1}
    caption = "Parallax factor coverage for an r=20 star (0 is bad, 0.5-1 is good). "
    caption += "One expects the parallax factor coverage to vary because stars on the ecliptic "
    caption += "can be observed when they have no parallax offset while stars at the pole are always "
    caption += "offset by the full parallax offset."""
    displayDict = {'group': parallaxGroup, 'subgroup': 'Parallax', 'order': order,
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1
    metric = metrics.ParallaxCoverageMetric(metricName='Parallax Coverage 24',
                                            rmag=24, seeingCol=geomSeeingCol)
    plotDict = {"xMin":     0, "xMax":     1, "binsize": 0.01,
                "colorMin": 0, "colorMax": 1}
    caption = "Parallax factor coverage for an r=24 star (0 is bad, 0.5-1 is good). "
    caption += "One expects the parallax factor coverage to vary because stars on the ecliptic "
    caption += "can be observed when they have no parallax offset while stars at the pole are always "
    caption += "offset by the full parallax offset."""
    displayDict = {'group': parallaxGroup, 'subgroup': 'Parallax', 'order': order,
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1
    metric = metrics.ParallaxDcrDegenMetric(metricName='Parallax-DCR degeneracy 20', rmag=20,
                                            seeingCol=geomSeeingCol)
    plotDict = {"xMin":     -0.2, "xMax":     0.2, "binsize": 0.01,
                "colorMin": -0.2, "colorMax": 0.2}
    caption = 'Correlation between parallax offset magnitude and hour angle an r=20 star.'
    caption += ' (0 is good, near -1 or 1 is bad).'
    displayDict = {'group': parallaxGroup, 'subgroup': 'Parallax', 'order': order,
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1
    metric = metrics.ParallaxDcrDegenMetric(metricName='Parallax-DCR degeneracy 24', rmag=24,
                                            seeingCol=geomSeeingCol)
    plotDict = {"xMin":     -0.2, "xMax":     0.2, "binsize": 0.01,
                "colorMin": -0.2, "colorMax": 0.2}
    caption = 'Correlation between parallax offset magnitude and hour angle an r=24 star.'
    caption += ' (0 is good, near -1 or 1 is bad).'
    displayDict = {'group': parallaxGroup, 'subgroup': 'Parallax', 'order': order,
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1

    metric = metrics.ProperMotionMetric(metricName='Proper Motion 20', rmag=20,
                                        seeingCol=geomSeeingCol)

    summaryStats = allStats
    plotDict = {'xMin':     0, 'xMax':     2, "binsize": 0.01,
                "colorMin": 0, "colorMax": 2}
    displayDict = {'group': properMotionGroup, 'subgroup': 'Proper Motion', 'order': order,
                   'caption': 'Proper Motion precision at r=20.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1
    metric = metrics.ProperMotionMetric(rmag=24, metricName='Proper Motion 24',
                                        seeingCol=geomSeeingCol)
    summaryStats = allStats
    plotDict = {'xMin':     0, 'xMax':     4, "binsize": 0.05,
                "colorMin": 0, "colorMax": 4}
    displayDict = {'group': properMotionGroup, 'subgroup': 'Proper Motion', 'order': order,
                   'caption': 'Proper Motion precision at r=24.'}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1
    metric = metrics.ProperMotionMetric(rmag=24, normalize=True, metricName='Proper Motion Normed',
                                        seeingCol=geomSeeingCol)
    plotDict = {'xMin':     0, 'xMax':     1, "binsize": 0.01,
                "colorMin": 0, "colorMax": 1}
    caption = 'Normalized proper motion at r=24. '
    caption += '(normalized to optimum observation cadence - start/end. 1=optimal).'
    displayDict = {'group': properMotionGroup, 'subgroup': 'Proper Motion', 'order': order,
                   'caption': caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)
    order += 1


    # And add a histogram of the time between quick revisits.
    order = 0

    binMin = 0
    binMax = 120.
    binsize = 3.
    bins_metric = np.arange(binMin / 60.0 / 24.0, (binMax + binsize) / 60. / 24., binsize / 60. / 24.)
    bins_plot = bins_metric * 24.0 * 60.0
    m1 = metrics.TgapsMetric(bins=bins_metric, metricName='dT visits')
    plotDict = {'bins': bins_plot, 'xlabel': 'dT (minutes)'}
    caption = ('Histogram of the time between consecutive revisits (<%.1f minutes), over entire sky.'
               % (binMax))
    displayDict = {'group': intraNightGroup, 'subgroup': 'Rapid Revisit', 'order': order,
                   'caption': caption}
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
    plotFunc = plots.SummaryHistogram()
    bundle = metricBundles.MetricBundle(m1, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, runName=runName,
                                        metadata=metadata, plotFuncs=[plotFunc])
    bundleList.append(bundle)
    order += 1

    # Loneliness Metric
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
    sqlconstraint = wfdWhere
    metadata = "WFD"

    cutoffMins = 60
    metric = LonelinessMetric(dTCutoff=cutoffMins*60)
    plotDict = {"xMin": 0, "xMax": 1, "binsize": 0.01,
                "colorMin": 0, "colorMax": 1,
                "xlabel": "Loneliness"}
    caption = "Fraction of visits that have no time-adjacent visit within " + \
              "{} minutes.".format(cutoffMins)
    displayDict = {"group": intraNightGroup, "subgroup": "Loneliness", "order": order,
                   "caption": caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)

    order += 1

    # Intranight Color Metric
    metric = IntranightColorMetric()
    plotDict = {"xMin": 0, "xMax": 1, "binsize": 0.01,
                "colorMin": 0, "colorMax": 1,
                "xlabel": "intra-night color fraction"}
    caption = "Considering only nights when this point was visited at least twice," + \
              " this metric returns the fraction of nights when a color was obtained" + \
              " at this point (i.e. visits in 2+ filters)"
    displayDict = {"group": intraNightGroup, "subgroup": "Intra-Night Color",
                   "order": order, "caption": caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName, metadata=metadata)
    bundleList.append(bundle)

    order += 1

    metric = SNLotsMetric()
    sqlconstraint = wfdWhere + ' and (filter="r" or filter="g" or filter="i" or filter="z")'
    plotDict = {"xMin": 0, "xMax": 0.4, "binsize": 0.01,
                "colorMin": 0, "colorMax": 0.4,
                "xlabel": "Fraction SN Ia well-sampled"}
    caption = "Fraction of z=0.5 type Ia SN that are detected in at least 3 filters" + \
              " before peak and twice in all filters post-peak (considering only" + \
              " the griz filters)"
    displayDict = {"group": snGroup, "subgroup": "SNLots", "order": order,
                   "caption": caption}
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, plotDict=plotDict,
                                        displayDict=displayDict, summaryMetrics=summaryStats,
                                        runName=runName)
    bundleList.append(bundle)
    order += 1
 
    # Median inter-night gap (each and all filters)
    slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
    metric = metrics.InterNightGapsMetric(metricName='Median Inter-Night Gap')
    sqls = ['filter = "%s"' % f for f in filters]
    orders = [filtorder[f] for f in filters]
    orders.append(0)
    sqls.append('')
    for sql, order in zip(sqls, orders):
        displayDict = {'group': interNightGroup, 'subgroup': 'Median Gap',
                       'caption': 'Median gap between days',
                       'order': order}
        if order == 0:
            xMin = 0
            xMax = 6
            binsize = 0.05
        elif order == filtorder["u"] or order == filtorder["g"]:
            xMin = 0
            xMax = 90
            binsize = 1
        else:
            xMin = 0
            xMax = 50
            binsize = 1

        plotDict = {"xMin": xMin, "xMax": xMax, "binsize": binsize,
                    "colorMin": xMin, "colorMax": xMax}
        bundle = metricBundles.MetricBundle(metric, slicer, sql, plotDict=plotDict,
                                            displayDict=displayDict, runName=runName)
        bundleList.append(bundle)

    # NEO XY plots
    slicer = slicers.UniSlicer()
    metric = metrics.PassMetric(metricName='NEODistances')
    stacker = stackers.NEODistStacker()
    stacker2 = stackers.EclipticStacker()
    for f in filters:
        plotFunc = plots.NeoDistancePlotter(eclipMax=10., eclipMin=-10.)
        caption = 'Observations within 10 degrees of the ecliptic. Distance an H=22 NEO would be detected'
        displayDict = {'group': neoGroup, 'subgroup': 'xy', 'order': filtorder[f],
                       'caption': caption}
        plotDict = {}
        sqlconstraint = 'filter = "%s"' % (f)
        bundle = metricBundles.MetricBundle(metric, slicer,
                                            sqlconstraint, displayDict=displayDict,
                                            stackerList=[stacker, stacker2],
                                            plotDict=plotDict, runName=runName,
                                            plotFuncs=[plotFunc])
        bundleList.append(bundle)




    return metricBundles.makeBundlesDictFromList(bundleList)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Python script to run MAF cadence performance metrics")
    parser.add_argument('dbFile', type=str, default=None,
            help="full file path to the opsim sqlite file")
    parser.add_argument("--outDir", type=str, default='./Out',
            help='Output directory for MAF outputs.')
    parser.add_argument("--runName", type=str, help="Name of the run")

    parser.set_defaults()
    args, _ = parser.parse_known_args()

    resultsDb = db.ResultsDb(outDir=args.outDir)
    opsdb = db.OpsimDatabaseV4(args.dbFile)
    bundleDict = makeBundleList(args.dbFile, runName=args.runName)
    group = metricBundles.MetricBundleGroup(bundleDict, opsdb, outDir=args.outDir,
                                            resultsDb=resultsDb)
    
    group.runAll()
    group.plotAll()

    opsdb.close()
