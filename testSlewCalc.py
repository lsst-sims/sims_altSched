from __future__ import division

import Telescope
import sqlite3
from matplotlib import pyplot as plt
import numpy as np

conn = sqlite3.connect("../opsim_runs/minion_1012_sqlite.db")
c = conn.cursor()

mySlewTimes = []
theirSlewTimes = []

prevAltaz = np.zeros(2)
prevNight = -1
prevFilter = 'r'

i = 0
for alt, az, filter, night, slewTime in c.execute("select altitude, azimuth, filter, night, slewTime from Summary group by expMJD"):
    altaz = np.array([alt, az])
    mySlewTime = Telescope.calcSlewTime(prevAltaz, np.array([alt, az]))

    if filter != prevFilter:
        mySlewTime = max(120, mySlewTime)

    if prevNight == night:
        mySlewTimes.append(mySlewTime)
        theirSlewTimes.append(slewTime)

    prevAltaz = altaz
    prevFilter = filter
    prevNight = night
    i += 1
    if i > 1000:
        break

mySlewTimes = np.array(mySlewTimes)
theirSlewTimes = np.array(theirSlewTimes)

print "my mean", np.mean(mySlewTimes), ", median", np.median(mySlewTimes), ", std", np.std(mySlewTimes)
print "their mean", np.mean(theirSlewTimes), ", median", np.median(theirSlewTimes), ", std", np.std(theirSlewTimes)

bins = [i/2 for i in range(260)]
plt.title("Histogram of slew times as calculated by opsim")
plt.xlabel("Slew time (sec)")
plt.hist(theirSlewTimes, bins, normed=1)

plt.figure()
plt.title("Cumulative distribution of slew times as calculated by opsim")
plt.xlabel("Slew time (sec)")
plt.hist(theirSlewTimes, 500, histtype="step", cumulative=True)

plt.figure()
plt.title("Cumulative sum of opsim slew times over slew rank")
sortedTimes = np.sort(theirSlewTimes)
cum = np.cumsum(sortedTimes)
plt.plot(np.arange(len(cum)), cum)

plt.figure()
plt.title("Cumulative sum of opsim slew times over slew rank (w/o filter change slews)")
sortedTimes = sortedTimes[sortedTimes != 120]
cum = np.cumsum(sortedTimes)
plt.plot(np.arange(len(cum)), cum)
print "excluding filters"
print "total cumulative slew time: ", cum[-1]
print "rank @ half total cum / # slews", np.searchsorted(cum, cum[-1]/2) / len(cum)

# show agreement between my calculated slew times and their calculated slew times
plt.figure()
plt.title("Agreement between my and their slew calculations")
plt.xlabel("Their slew time (sec)")
plt.ylabel("My slew time (sec)")
plt.scatter(theirSlewTimes, mySlewTimes)

minTime = min(mySlewTimes)
maxTime = max(mySlewTimes)
plt.plot([0, maxTime], [0, maxTime])
plt.xlim(0, maxTime)
plt.ylim(0, maxTime)


plt.show()

