from __future__ import division, print_function
import numpy as np
from astropy import wcs
import glob
import os
import utils

class Tiling:
    def getTiling(self, density):
        # density is in pointings / steradian
        raise NotImplementedError("Subclass has not implemented getTiling")

class RandomTiling(Tiling):
    def getTiling(self, density):
        numPointings = int(4 * np.pi * density)
        thetas = np.arcsin(np.random.random(numPointings)*2 - 1)
        phis = np.random.random(numPointings)*2*np.pi

        pointings = np.vstack([phis, thetas]).transpose()
        return pointings


class TwoRaftTiling(Tiling):
    def getTiling(self, density):
        """ Tiles the sky with a 2-raft-overlap pattern

        Or rather, the tiling has 2 rafts of overlap in flat space. Turns
        out that after projecting onto the sphere, this isn't really the
        case any more. Fortunately, since we do large-scale full-tiling
        dithering, the exact tiling doesn't matter all that much.
        """
        # to generate the pointings, start with the pointing pattern in flat
        # space, then map the pattern onto a quadrilateralized spherical cube

        # the coordinates in flat space are (u,v) where u and v represent
        # the number of rafts of displacement in the horizontal and vertical
        # directions, respectively, from the bottom left corner (0,0)

        # in these coordinates, focal plane is a square with side length
        # equal to 5 (since the focal plane is 5x5 rafts)

        # although the raft and fov widths are not hard-coded, the number of 
        # rafts is hard-coded. After all, this algorithm completely relies on
        # the fact that the focal plane is a 5x5 square with 1x1 squares
        # removed from the corners

        # I'm guessing there's a better way to do this, but here goes...

        #raftWidth = telescope.raftWidth * 0.95
        # TODO this is hacky but it gives about the right number
        # of pointings given the input density (the reason it's coded
        # in terms of raft widths is because that's how I originally
        # wrote it, but I don't anticipate this tiling being super
        # useful anyway, so haven't rewritten it in terms of density)
        raftWidth = np.sqrt(1 / (27 * density))

        cubeSideLength = np.arccos(1/3) # in radians
        numRaftsPerCubeSide = int(cubeSideLength / raftWidth)

        minU = 0
        minV = 0

        # we need to tile enough to be sure that we can completely cover
        # the unfolded cube with the pattern
        #maxU = int(2 * np.pi / raftWidth)
        #maxV = int(np.pi / raftWidth)

        # if you unfold the cube like this:
        #  _
        # | |_____
        # |  _____|
        # |_|
        # then we need pattern a rectangle of size numRaftsPerCubeSide * (4,3)
        # I'm using (5,4) instead to avoid the edges of the pattern
        maxU = 5 * numRaftsPerCubeSide 
        maxV = 4 * numRaftsPerCubeSide

        u = minU
        v = minV
        pointings = []

        # this while loop populates the pointings array with (u,v) tuples 
        # in the proper pointing pattern. The pattern fills a rectangle
        # stretching from (minU, minV) to (maxU, maxV)
        while True:
            prevBottomV = v
            prevBottomU = u
            while True:
                # scan up and to the left, moving left 1 raft and up 4 rafts 
                # each time

                pointings.append((u,v))

                u -= 1
                v += 4 

                if u < minU or v > maxV:
                    # we're done with this column
                    break
                    
            if prevBottomV >= 1: 
                # start the next column shifted 5 rafts right and one raft down
                u = prevBottomU + 5 
                v = prevBottomV - 1
            else:
                # start the next column shifted 4 rafts right and 3 rafts up
                u = prevBottomU + 4
                v = prevBottomV + 3

            # now check if we've reached the bottom right corner
            # if so, we need to start the next column farther up
            if u > maxU:
                v += 4 * (u - maxU)
                u = maxU

                # if v is now larger than maxV, that means we're done
                if v > maxV:
                    break

        pointings = np.array(pointings).astype("float64")

        
        # astropy lays out the cube like this
        # _
        #|0|_____
        #|1 4_3_2|
        #|5|
        # 
        # where (0,0) is (annoyingly) in the center of face 1
        # and where faces 4, 3, and 2 are repeated to the left of
        # face 1
        #     

        # convert to "degrees." Each cube face has a side length of 90 degrees 
        pointings *= 90 / numRaftsPerCubeSide

        # shift the origin so the pattern covers the entire unfolded cube
        # (and doesn't double cover anything)
        pointings[:,0] -= 45
        pointings[:,1] -= 135 

        # convert to RA/DEC
        w = wcs.WCS(naxis=2)
        # QSC is quadrilateralized spherical cube
        w.wcs.ctype = ["RA---QSC", "DEC--QSC"]
        # the 2nd argument is the "origin" of the coordinate system
        # "1" represents (1,1) and I think it conveys that the QSC 
        # point (0,0) corresponds to RA/DEC of (0,0)? TODO
        pointings = w.wcs_pix2world(pointings, 0) 
        
        # nan indicates the pointing was not on any cube face
        pointings = pointings[~np.isnan(pointings[:,0])]

        # and convert back to radians
        pointings = np.radians(pointings)

        return pointings


class ThomsonTiling(Tiling):
    def __init__(self):
        # check which Thomson problem solution files are available
        solFiles = glob.glob(os.path.join("thomsonSols","*.txt"))
        solFiles = [os.path.basename(solFile) for solFile in solFiles]
        self.availN = np.array([int(solFile[:-4]) for solFile in solFiles])

    def getTiling(self, density):
        numPointings = int(4 * np.pi * density)
        N = self.availN[np.argmin(np.abs(self.availN - numPointings))]
        xyz = np.zeros((N,3))
        with open(os.path.join("thomsonSols", str(N) + ".txt"), "r") as solFile:
            for i, line in enumerate(solFile):
                xyz[i,:] = list(map(float, line.split()))
        phiTheta = utils.cartesian2Spherical(xyz[:,0], xyz[:,1], xyz[:,2])
        return phiTheta
