from __future__ import division

import matplotlib
# use this backend because it supports the Yield function
# thus allowing the user to interact with the figure
# while we keep computing the next Visit
matplotlib.use("wxAgg")
import wx

from mayavi import mlab
import numpy as np
import time

class GraphicalMonitor3D:
    def __init__(self, context, res):
        self.context = context
        self.res = res

    def init(self):
        # TODO parametrize

        # the step size complex to make the end of the range inclusive
        phi, theta = np.mgrid[0:2*np.pi:2*self.res*1j, -np.pi/2:np.pi/2:self.res*1j]

        x = np.cos(theta) * np.cos(phi)
        y = np.cos(theta) * np.sin(phi)
        z = np.sin(theta)

        self.skyMap = np.zeros_like(x)
        self.i = 0

        self.skyPlt = mlab.mesh(x, y, z, scalars = np.zeros_like(x))
        return


    def addVisit(self, visit):
        # update self.skyMap with the new informations
        # round visit coords to nearest sky pixel
        phiId = int(visit.ra / (2 * np.pi) * self.res)
        thetaId = int((visit.dec + np.pi / 2)/ np.pi * self.res)
        # now find all pixels contained in the shape by projecting
        # each candidate pixel onto the plane tangent to the sphere
        # at (visit.ra, visit.dec)

        self.skyMap[phiId, thetaId] += 1
        wx.Yield()

    def updateDisplay(self):
        self.skyPlt.mlab_source.set(scalars = self.skyMap)
        wx.Yield()
