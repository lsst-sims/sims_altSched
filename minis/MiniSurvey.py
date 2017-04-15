from  __future__ import division
import numpy as  np
import Telescope
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from RotationGenerator import RotationGenerator
from Visit import VisitPair
import Config
import Utils
from astropy import wcs

class MiniSurvey:
    NORTH = 0
    SOUTH = 1

    # the very first mini should have rotation -pi/2
    prevRot = -1
    prevPrevRot = -1

    rotationGenerator = RotationGenerator()

    def __init__(self, pointings, rotation):
        if not isinstance(self, MiniSurvey):
            raise NotImplementedError("call MiniSurvey.newMiniSurvey() to \
                                       get a new MiniSurvey instance")

        # the VisitPair constructor relies on self.rotation having been
        # set already. TODO do something that's not this
        self.rotation = rotation

        self.visitPairs = [VisitPair(self, pointing[0], pointing[1])
                           for pointing in pointings]


        self.numIncompleteVisitPairs = len(self.visitPairs)

    @classmethod
    def setLatestRotation(cls, rotation, direction):
        # called when resuming the survey from a checkpoint
        cls.rotationGenerator = RotationGenerator(rotation, direction)

    @classmethod
    def newMiniSurvey(cls, minDec, maxDec, minRa, maxRa):
        rotation = next(cls.rotationGenerator.rotations())
        #allPointings = cls._generateRandomPointings()
        allPointings = cls._realGeneratePointings()

        # now take the pointings and rotate them around the z and x
        # axes so they are randomly dithered on the sphere
        zDither = np.random.randn()
        xDither = np.random.randn()
        
        # rotating around z is easy
        allPointings[:,0] += zDither

        # rotating around x is harder
        sin = np.sin
        cos = np.cos
        phi = allPointings[:,0]
        theta = allPointings[:,1]
        newPhi = np.arctan2(-sin(xDither)*sin(theta) + cos(xDither)*cos(theta)*sin(phi),
                            cos(theta)*cos(phi))
        newPhi %= 2 * np.pi
        newTheta = np.arcsin(cos(xDither)*sin(theta) + sin(xDither)*cos(theta)*sin(phi))
        allPointings = np.vstack([newPhi, newTheta]).T

        #print "RA: (", raStart, ",", raEnd, ")"
        #print "DEC: (", decStart, ",", decEnd, ")"

        # this mini survey should only include pointings between ra/dec min/max
        pointings = allPointings[np.where(
            Utils.areRasInRange(allPointings[:,0], (minRa, maxRa)) & 
            (minDec < allPointings[:,1]) & (allPointings[:,1] < maxDec)
        )]

        if len(pointings) > 1000:
            # TODO debug: there shouldn't ever be more than 1000 pointings
            # in a mini survey, so this will show you what the pointings
            # are if this happens
            fig = plt.figure()
            ax = Axes3D(fig)
            pointings = np.array(pointings)
            phis = pointings[:,0]
            thetas = pointings[:,1]
            ax.scatter(np.cos(phis)*np.cos(thetas),
                       np.sin(phis)*np.cos(thetas),
                       np.sin(thetas))

            plt.show()

        newMini = cls(pointings, rotation)
        return newMini
        

    @classmethod
    def _generateRandomPointings(self):
        numPointings = int(4 * np.pi / (Telescope.fovWidth**2))
        thetas = np.arcsin(np.random.random(numPointings)*2 - 1)
        phis = np.random.random(numPointings)*2*np.pi

        pointings = np.vstack([thetas, phis]).transpose()
        return pointings


    @classmethod
    def _realGeneratePointings(self):
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

        cubeSideLength = np.arccos(1/3) # in radians
        numRaftsPerCubeSide = int(cubeSideLength / Telescope.raftWidth)

        minU = 0
        minV = 0

        # we need to tile enough to be sure that we can completely cover
        # the unfolded cube with the pattern
        #maxU = int(2 * np.pi / Telescope.raftWidth)
        #maxV = int(np.pi / Telescope.raftWidth)

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

        # offset by the amount in the arguments
        # TODO no longer doing this
        #pointings[:,0] += raOffset
        #pointings[:,1] += decOffset

        #plt.scatter(pointings[:,0],pointings[:,1])
        #plt.show()
        #exit()


        #print "RA: (",min(pointings[:,0]),",",max(pointings[:,0]),")"
        #print "DEC:(",min(pointings[:,1]),",",max(pointings[:,1]),")"

        return pointings

        """
        cubeFaceLowerLeftCorners = np.array([[1,2], [0,1], [1,1], 
                                             [2,1], [3,1], [1,0]])
        cubeFaceUpperRightCorners = cubeFaceLowerLeftCorners + 1

        cubeFaceLowerLeftCorners *= numRaftsPerCubeSide
        cubeFaceUpperRightCorners *= numRaftsPerCubeSide

        # add 6 rafts to u and v avoid the edges
        cubeFaceLowerLeftCorners += 6
        cubeFaceUpperRightCorners += 6

        thetas = []
        phis = []
        for pointing in pointings:
            (u, v) = pointing
            face = -1
            for i in range(6):
                if (np.all(pointing > cubeFaceLowerLeftCorners[i]) and
                    np.all(pointing < cubeFaceUpperRightCorners[i])):
                    face = i

            if 1 <= face <= 4:
                mu = np.arctan(v / np.sqrt(numRaftsPerCubeSide**2+u**2))
                # use (face-2) so face2 gets the point (mu, nu) = (0,0)
                # TODO center it to minimize distortion in the actual minisurvey
                # area and maybe dither randomly
                nu = np.arctan(u / numRaftsPerCubeSide) + (face - 2) * np.pi/2

            elif face == 0:
                mu = np.arctan(numRaftsPerCubeSide / np.sqrt(u**2+v**2))
                nu = np.arctan2(u, -v)

            elif face == 5:
                mu = -np.arctan(numRaftsPerCubeSide / np.sqrt(u**2+v**2))
                nu = np.arctan2(-u, v)

            else:
                # this pointing doesn't fall within the cube
                continue

            if face != 5:
                continue

            t = np.pi / 12 * np.tan(mu)
            theta = np.arctan(np.sin(t) / (np.cos(t) - (1 / np.sqrt(2))))
            cosphi = ( (1 - np.cos(mu)**2 * np.tan(nu)**2) *
                       (1 - np.cos(np.arctan(1/np.cos(theta)))) )
            phi = np.arccos(cosphi)
            
            thetas.append(theta)
            phis.append(phi)

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(np.cos(thetas)*np.cos(phis), 
                   np.sin(thetas)*np.cos(phis), 
                   np.sin(phis))
        plt.show()
        
        pointings *= Telescope.raftWidth
        
        ra += raOffset
        dec += decOffset

        visitPairs = []
        for i in range(len(ra)):
            visitPairs.append(VisitPair(self, ra[i], dec[i]))

        return pointings
        """
