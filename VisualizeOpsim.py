from __future__ import division

import sqlite3
import Config
from graphics.GraphicalMonitor import GraphicalMonitor
import astropy.time
from minis.Visit import Visit
from Simulator import Simulator
import AstronomicalSky
from SummaryPlots import SummaryPlots
from minis.SkyMap import SkyMap

from Telescope import Telescope

import time
from matplotlib import pyplot as plt
import sys
import numpy as np

showDisp = True
clearDisplayNightly = True

class VisualizeOpsim:

    def __init__(self, imScale=2):
        conn = sqlite3.connect("../opsim_runs/minion_1012_sqlite.db")
        self.c = conn.cursor()
        self.imScale = imScale
        self.curTime = Config.surveyStartTime


    def main(self):
        tel = Telescope()
        if showDisp:
            skyMap = SkyMap(telescope = tel, resScale=2)
            display = GraphicalMonitor(skyMap=skyMap)

        i = 0
        prevI = i
        isNightYoung = True
        prevNightNum = -1

        alts = []
        azes = []
        obsDecs = []

        for row in self.c.execute("select expMJD, night, ditheredRA, ditheredDEC, visitTime from ObsHistory"):
            mjd = row[0]
            nightNum = row[1]
            ra = row[2]
            dec = row[3]
            expTime = row[4]

            t = astropy.time.Time(mjd, format="mjd")
            self.curTime = t.unix

            radec = np.array([[ra, dec]])
            altaz = AstronomicalSky.radec2altaz(radec, self.curTime)
            alts.append(altaz[0,0])
            azes.append(altaz[0,1])
            obsDecs.append(dec)

            # no slew time is over an hour
            isNightYoung = nightNum > prevNightNum
            if isNightYoung:
                print "Night:", nightNum, "\r",
                sys.stdout.flush()

            visit = Visit(None, ra, dec, 0, expTime)

            if showDisp:
                skyMap.addVisit(visit, self.curTime)

            perNight, deltaI = Simulator.getUpdateRate(i)

            if showDisp and ((perNight and isNightYoung) or 
                             (not perNight and i - prevI >= deltaI)):
                display.updateDisplay(skyMap, self.curTime)
                #display.saveFrame("images/opsim/%07d.png" % i)
                prevI = i

            if isNightYoung and showDisp and clearDisplayNightly:
                skyMap.clear()
            i += 1
            prevNightNum = nightNum
            if nightNum > 30:
                break

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
        plt.show()

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
    V = VisualizeOpsim(imScale=10)
    V.main()
