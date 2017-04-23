from __future__ import division

import sqlite3
from SummaryPlots import SummaryPlots

conn = sqlite3.connect("../opsim_runs/minion_1012_sqlite.db")
c = conn.cursor()
altazdec = c.execute("select altitude, azimuth, fieldDec from Summary group by expMJD limit 1000")

# this syntax "unzips" altazdec
alt, az, dec = map(list, zip(*list(altazdec)))

plotter = SummaryPlots(alts=alt, azes=az, decs=dec)
plotter.dAirmassContour()
plotter.dAirmassCum()
plotter.zenithAngleContour()
plotter.show()
