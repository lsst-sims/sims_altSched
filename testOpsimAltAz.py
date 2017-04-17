from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
import sqlite3
import Telescope


def showAirmassPlots(plotter, altazdec):

    # dAirmass = delta airmass = observed airmass - optimal airmass given dec
    # first, create a polar plot with r = dAirmass and theta = azimuth

    # limits and resolution for r and theta
    nDAirmass = 100
    minDAirmass = 0
    maxDAirmass = 0.5
    nAzes = 150

    r, theta = np.meshgrid(np.linspace(minDAirmass, maxDAirmass, num=nDAirmass),
                           np.linspace(0, 2*np.pi, num=nAzes))
    values = np.zeros(r.shape)
    
    # keep track of the airmasses and dAirmasses bserved at so we can make
    # histograms
    obsAirmasses = []
    dAirmasses = []

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
        
        # don't plot points out of range
        if dAirmass > maxDAirmass or dAirmass < minDAirmass:
            continue

        # rId and thetaId are indices within r and theta corresponding
        # to the dAirmass and azimuth values
        rId = (dAirmass-minDAirmass) / (maxDAirmass - minDAirmass) * nDAirmass
        rId = int(rId)
        thetaId  = int(az / (2*np.pi) * nAzes)
        values[thetaId, rId] += 1

    # first, make a histogram showing the cumulative distribution of dAirmasses
    plotter.figure()
    plotter.title("Cumulative distribution of observed airmass - optimal airmass")
    plotter.xlabel("Number of airmasses observed through - " + \
                   "sec(declination - telescope latitude)")
    plotter.hist(dAirmasses, 300, cumulative=True)

    # now make a histogram of the airmass observed through
    plotter.figure()
    plotter.title("Histogram of airmasses observed through")
    plotter.xlabel("airmass = csc(altitude)")
    plotter.hist(obsAirmasses, 300)

    # now show a contour plot in polar coordinates where r is dAirmass
    # and theta is azimuth
    fig, ax = plotter.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    ax.set_ylim(minDAirmass, maxDAirmass)
    p = ax.contourf(theta, r, np.log(values))
    plotter.colorbar(p, ax=ax)

if __name__ == "__main__":
    conn = sqlite3.connect("../opsim_runs/minion_1012_sqlite.db")
    c = conn.cursor()
    altazdec = c.execute("select altitude, azimuth, fieldDec from Summary group by expMJD")
    showAirmassPlots(plt, altazdec)
    plt.show()
