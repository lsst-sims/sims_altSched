from __future__ import division

from graphics.GraphicalMonitor import GraphicalMonitor
#from graphics.GraphicalMonitor3D import GraphicalMonitor3D

import numpy as np
from minis.Scheduler import Scheduler
import Config
import AstronomicalSky
import Telescope
import time

from matplotlib import pyplot as plt

showDisp = True
dispUpdateFreq = 1
saveMovie = False

# TODO plot azimuth vs time
# change filters every revisit set (revisits 1 hour, change every 2 hours)
# add in settle time

class Simulator:
    def __init__(self):
        self.startTime = Config.surveyStartTime
        self.curTime = self.startTime

    def run(self):
        sched = Scheduler(context=self)

        # imScale is the resolution (and therefore the speed) that we want
        if showDisp:
            display = GraphicalMonitor(context=self, imScale=4, numVisitsDisplayed=None)
            #display = GraphicalMonitor3D(context = self, res = 180)
            #display.init()
       
        nightNum = 0
        isNightYoung = True

        slewTimes = []
        revisitTimes = []

        i = 0
        for visit in sched.schedule():
            # if visit is None, that means the scheduler ran out of places
            # to point for the night
            if visit is not None:
                radec = np.array([[visit.ra, visit.dec]])
                altaz = AstronomicalSky.radec2altaz(radec, self.time())[0]
                if not isNightYoung:
                    slewTime = Telescope.calcSlewTime(prevAltaz, altaz)

                if showDisp:
                    display.addVisit(visit)      

            if showDisp:
                if i % dispUpdateFreq == 0:
                    display.updateDisplay()
                    #plt.savefig("images/pernight/%04d.png" % nightNum)
                    if saveMovie:
                        display.saveFrame("images/pygame/%05d.png" % (i / dispUpdateFreq))

            if i % 10000 == 0:
                print i

            
            if visit is None:
                self.curTime += 30
            else:
                # add the exposure time of this visit to the current time
                expTime = AstronomicalSky.getExpTime(visit.ra, visit.dec)
                self.curTime += expTime
                if not isNightYoung:
                    self.curTime += slewTime
                    slewTimes.append(slewTime)
                sched.notifyVisitComplete(visit, self.time())
                visitPair = visit.visitPair
                if visitPair.visit1.isComplete and visitPair.visit2.isComplete:
                    revisitTimes.append(np.abs(visitPair.visit1.timeOfCompletion - \
                                               visitPair.visit2.timeOfCompletion))
                prevAltaz = altaz

            if self.curTime < AstronomicalSky.nightEnd(nightNum):
                isNightYoung = False
            else:
                self.curTime = AstronomicalSky.nightStart(nightNum + 1)
                nightNum += 1
                isNightYoung = True

            i += 1
            #if showDisp and saveVid:
            if i > 30000:
                print "avg slew time", np.mean(slewTimes), "seconds"
                plt.hist(slewTimes, bins = np.arange(min(slewTimes), max(slewTimes), 0.5))
                plt.xlabel("Slew Time (secs)")
                plt.ylabel("Number of slews")
                plt.title("Histogram of Slew Times")

                revisitTimesMins = np.array(revisitTimes) / 60
                print "avg revisit time", np.mean(revisitTimesMins), "minutes"
                plt.figure()
                plt.hist(revisitTimesMins, 
                         bins = np.arange(0, 3*np.median(revisitTimesMins), 1))
                plt.xlabel("Revisit Time (mins)")
                plt.ylabel("Number of Revisits")
                plt.title("Histogram of Revisit Times")

                plt.show()
                return

    def time(self):
        return self.curTime

if __name__ == "__main__":
    sim = Simulator()
    sim.run()
