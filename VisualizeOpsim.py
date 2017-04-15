from __future__ import division

import sqlite3
import Config
from graphics.GraphicalMonitor import GraphicalMonitor
import astropy.time
from minis.Visit import Visit
from Simulator import Simulator

import time

class VisualizeOpsim:

    def __init__(self, imScale=2):
        conn = sqlite3.connect("../opsim_runs/minion_1012_sqlite.db")
        self.c = conn.cursor()
        self.imScale = imScale
        self.curTime = Config.surveyStartTime


    def main(self):
        display = GraphicalMonitor(context=self, imScale=self.imScale)

        i = 0
        prevI = i
        isNightYoung = True
        prevNightNum = -1

        for row in self.c.execute("select expMJD, night, ditheredRA, ditheredDEC, visitTime from ObsHistory"):
            mjd = row[0]
            nightNum = row[1]
            ra = row[2]
            dec = row[3]
            expTime = row[4]

            t = astropy.time.Time(mjd, format="mjd")
            self.curTime = t.unix

            # no slew time is over an hour
            isNightYoung = nightNum > prevNightNum

            visit = Visit(None, ra, dec, 0, expTime)

            display.addVisit(visit)

            perNight, deltaI = Simulator.getUpdateRate(i)

            if (perNight and isNightYoung) or (not perNight and i - prevI >= deltaI):
                display.updateDisplay()
                display.saveFrame("images/opsim/%07d.png" % i)
                prevI = i

            i += 1
            prevNightNum = nightNum
            

    def time(self):
        return self.curTime

if __name__ == "__main__":
    V = VisualizeOpsim(imScale=5.5)
    V.main()
