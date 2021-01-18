import numpy
import math
import random
import MHDCoordSys
from scipy.spatial.transform import Rotation as R
import MHDRectPM

class MHDMagAssby:
    def __init__(self, csIn):
        self.cs = csIn
        self.magnets = []

    def addMagnet(self, magIn):
        self.magnets.append(magIn)

    def getBForI(self, I_MagAssby):
        # I_MagAssem is a vector in MagAssem's CS pointing from MagAssem's center to magnetometer
        bSum_MagAssby = numpy.array([0.0, 0.0, 0.0], numpy.double)
        for mag in self.magnets:
            # lN vector points from MagAssem's center to Nth mag's center
            lN_MagAssby = numpy.array([mag.cs.x, mag.cs.y, mag.cs.z], numpy.double)
            # Compute vector from magnet to magnetomoter
            mN_MagAssby = I_MagAssby - lN_MagAssby
            # Rotate m1 into mag coords
            rotMag = R.from_euler(mag.cs.rot, [mag.cs.thetaX, mag.cs.thetaY, mag.cs.thetaZ], degrees=True)
            rInvMag = rotMag.inv()
            mN_Mag = rInvMag.apply(mN_MagAssby)
            # compute magnetic field at point in magnet's CS
            bN_Mag = mag.getB(mN_Mag)
            # rotate to magAssby's CS
            bN_MagAssby = rotMag.apply(bN_Mag)
            # add to sum
            bSum_MagAssby+=bN_MagAssby

        return bSum_MagAssby

