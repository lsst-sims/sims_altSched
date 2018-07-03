from __future__ import division

import sqlite3
import config
from graphics.GraphicalMonitor import GraphicalMonitor
import astropy.time
from Visit import Visit
from Visit import PROP_WFD
from Simulator import Simulator
from lsst.sims.speedObservatory import sky
from SummaryPlots import SummaryPlots
from SkyMap import SkyMap

from lsst.sims.speedObservatory import Telescope

import time
from matplotlib import pyplot as plt
import sys
import numpy as np

showDisp = True
clearDisplayNightly = True

class VisualizeOpsim:

    def __init__(self, dbFilename, imScale=4):
        conn = sqlite3.connect(dbFilename)
        self.c = conn.cursor()
        self.imScale = imScale
        self.curTime = config.surveyStartTime


    def main(self):
        tel = Telescope()
        if showDisp:
            skyMap = SkyMap(telescope = tel, resScale=self.imScale)
            display = GraphicalMonitor(skyMap=skyMap, mode="filters")

        i = 0
        prevI = i
        isNightYoung = True
        prevNightNum = -1

        alts = []
        azes = []
        obsDecs = []

        for row in self.c.execute("select observationStartMJD, night, fieldRA, fieldDec, visitExposureTime, filter, slewTime from SummaryAllProps"):
            mjd = row[0]
            nightNum = row[1]
            ra = np.radians(row[2])
            dec = np.radians(row[3])
            expTime = row[4]
            filter = row[5]
            slewTime = row[6]

            t = astropy.time.Time(mjd, format="mjd")
            self.curTime = t.unix

            radec = np.array([[ra, dec]])
            alt, az = sky.radec2altaz(ra, dec, self.curTime)
            alts.append(alt)
            azes.append(az)
            obsDecs.append(dec)

            # no slew time is over an hour
            isNightYoung = nightNum > prevNightNum
            if isNightYoung:
                print("Night:", nightNum, end="\r")
                sys.stdout.flush()

            visit = Visit(PROP_WFD, None, ra, dec, 0, expTime, filter=filter)

            if showDisp:
                skyMap.addVisit(visit, slewTime, self.curTime)

            perNight, deltaI = False, 1 #Simulator.getUpdateRate(i / 1000, i)

            if showDisp and ((perNight and isNightYoung) or 
                             (not perNight and i - prevI >= deltaI)):
                display.updateDisplay(skyMap, self.curTime)
                #display.saveFrame("images/opsim/%07d.png" % i)
                prevI = i

            if isNightYoung and showDisp and clearDisplayNightly:
                skyMap.clear()
            i += 1
            prevNightNum = nightNum
            if nightNum > 365 - 1:
                break

        return

        # code below here makes plots, but may no longer work
        revisitTimesMap = skyMap.getRevisitMap()
        allRevisitTimes = []
        for pix in revisitTimesMap.flatten():
            if isinstance(pix, list):
                allRevisitTimes += pix
        plt.figure("Revisit Time Histogram")
        plt.hist(np.array(allRevisitTimes)/3600, 300, cumulative=True, normed=True, range=[0,2])
        plt.title("Per-Pixel Revisit Times (opsim)")
        plt.xlabel("Time (hours)")
        plt.ylabel("Cumulative number of visits (normalized)")

        plt.figure("% Visits not accompanied by a revisit within 45 minutes")
        percentLonelyMap = skyMap.getLonelinessMap(cutoffMins=45)
        plt.title("Percent of visits with no revisit within 45 minutes")
        plt.imshow(percentLonelyMap)
        plt.clim(0,0.4)
        plt.colorbar()

        for percentile in [10, 50, 75, 90, 95, 99]:
            plt.figure(str(percentile) + "th percentile revisit time")
            percentileMap = skyMap.getPercentileMap(percentile)
            plt.title(str(percentile) + "th percentile revisit time (in days)")
            plt.imshow(percentileMap / 3600 / 24)
            plt.clim(0,7)
            plt.colorbar()


        plt.figure("avg revisit time")
        avgRevisitMap = skyMap.getAvgRevisitMap()
        plt.imshow(avgRevisitMap)
        plt.colorbar()
        plt.show(block=False)

        azes = np.array(azes) % (2*np.pi)

        plotter = SummaryPlots(alts=alts, azes=azes, decs=obsDecs)
        plotter.dAirmassCum()
        plotter.airmassHist()
        plotter.dAirmassContour()
        plotter.zenithAngleContour()
        plotter.show()
            

    def time(self):
        return self.curTime

if __name__ == "__main__":
    assert(len(sys.argv) == 2 and len(sys.argv[1]) > 0)
    dbFilename = sys.argv[1]
    V = VisualizeOpsim(dbFilename, imScale=3)
    V.main()
