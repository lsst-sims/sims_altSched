from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
import sqlite3
import Telescope


def showAirmassPlots(plotter, altazdec):

    # dAirmass = delta airmass = observed airmass - optimal airmass given dec
    # first, create a polar plot with r = dAirmass and theta = azimuth

    # limits and resolution for dAirmass
    nDAirmass = 50
    minDAirmass = 0
    maxDAirmass = 0.5

    # for a contour plot of zenith angle
    nZenith = 50
    minZenith = 0
    maxZenith = 90

    # limits for az are always 0 to 2pi
    # doesn't make sense to have resolution higher than ~3 degrees
    # since the FoV size will smear out the graph on that scale
    nAzes = 100

    # r, theta are polar coordinates in the plot
    DAr, theta = np.meshgrid(np.linspace(minDAirmass, maxDAirmass, num=nDAirmass),
                             np.linspace(0, 2*np.pi, num=nAzes))

    Zr, theta = np.meshgrid(np.linspace(minZenith, maxZenith, num=nZenith),
                            np.linspace(0, 2*np.pi, num=nAzes))
    dAirmassValues = np.zeros(DAr.shape)
    zenithValues = np.zeros(Zr.shape)
    
    # keep track of the airmasses and dAirmasses bserved at so we can make
    # histograms
    obsAirmasses = []
    dAirmasses = []
    zeniths = []

    # altazdec looks like [(alt1, az1, dec1), (alt2, az2, dec2), ...]
    for alt, az, dec in altazdec:
        # airmass number is csc of elevation or sec of zenith angle
        airmass = 1 / np.sin(alt)

        # the best place to observe a field is when it hits the meridian,
        # at which point you are observing through an airmass of
        # sec(dec - latitude)
        dAirmass = airmass - (1/np.cos(dec - Telescope.latitude))
        
        obsAirmasses.append(airmass)
        dAirmasses.append(dAirmass)
        zeniths.append(90 - np.degrees(alt))
        
        # don't plot points out of range
        if dAirmass > maxDAirmass or dAirmass < minDAirmass:
            continue

        # rId and thetaId are indices within r and theta corresponding
        # to the dAirmass/zenith and azimuth values
        DArId = (dAirmass-minDAirmass) / (maxDAirmass - minDAirmass) * nDAirmass

        # use this instead if you want to plot zenith angle instead of dAirmass
        ZrId = (90-np.degrees(alt) - minZenith) / (maxZenith - minZenith) * nZenith

        DArId = int(DArId)
        ZrId = int(ZrId)
        thetaId  = int(az / (2*np.pi) * nAzes)
        dAirmassValues[thetaId, DArId] += 1
        zenithValues[thetaId, ZrId] += 1

    # first, make a histogram showing the cumulative distribution of dAirmasses
    plotter.figure()
    plotter.title("Cumulative distribution of observed airmass - optimal airmass")
    plotter.xlabel("Number of airmasses observed through - " + \
                   "sec(declination - telescope latitude)")
    plotter.hist(dAirmasses, 300, range=(0,1), cumulative=True, normed=True)

    # now make a histogram of the airmass observed through
    plotter.figure()
    plotter.title("Histogram of airmasses observed through")
    plotter.xlabel("airmass = csc(altitude)")
    plotter.hist(obsAirmasses, 300)

    # now show contour plots in polar coordinates where r is dAirmass
    # or zenith angle and theta is azimuth
    fig, ax = plotter.subplots(subplot_kw=dict(projection='polar'))
    ax.set_title("log_10(nVisits)")
    ax.set_theta_zero_location("N")
    ax.set_ylim(minDAirmass, maxDAirmass)
    ax.set_xlabel("Azimuth (N=0deg, E=90deg)")
    ax.set_ylabel("Delta Airmass (obsAirmass - sec(obsDec - Tel.lat))")
    p = ax.contourf(theta, DAr, np.log10(dAirmassValues))
    plotter.colorbar(p, ax=ax)

    fix, ax = plotter.subplots(subplot_kw=dict(projection='polar'))
    ax.set_title("nVisits")
    ax.set_theta_zero_location("N")
    ax.set_ylim(minZenith, maxZenith)
    ax.set_xlabel("Azimuth (N=0deg, E=90deg)")
    ax.set_ylabel("Zenith angle (degrees)")
    p = ax.contourf(theta, Zr, zenithValues)
    plotter.colorbar(p, ax=ax)

if __name__ == "__main__":
    conn = sqlite3.connect("../opsim_runs/minion_1012_sqlite.db")
    c = conn.cursor()
    altazdec = c.execute("select altitude, azimuth, fieldDec from Summary group by expMJD limit 100000")
    showAirmassPlots(plt, altazdec)
    plt.show()
