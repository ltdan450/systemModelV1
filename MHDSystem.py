import MHDCoordSys
import MHDMagnetometer
#import MHDRectPM
#import MHDUniverse
import math
import numpy
from scipy.spatial.transform import Rotation as R

class MHDSystemA:

    def __init__(self,cs,universe):

        self.cs = cs
        self.universe = universe
        self.magnetometers=[]
        self.magnet = 0

        print ('need to setup self.magnetModel in MHDSystem')

    def addMagnet(self,magnet):
        self.magnet = magnet

    def addMagnetometer(self,mgtr):
        self.magnetometers.append(mgtr)

    def step(self,t):
        print('system t:%s'%t)

        oValid = False
        oIterations = 0
        oError = 10.0
        oErrorThresh = 0.001

        #while(oError>oErrorThresh):
        if 1 == 1:
            print('o iterations:%s'%oIterations)
            for mgtr in self.magnetometers:
                bMeasMgtr = mgtr.measure()      #what mgtr measures in its own coordinate system
                rMgtr = R.from_euler(mgtr.cs.rot, [mgtr.cs.thetaX, mgtr.cs.thetaY, mgtr.cs.thetaZ], degrees=True)
                bMeasSys = rMgtr.apply(bMeasMgtr)
                rMag = R.from_euler(self.magnet.cs.rot, [self.magnet.cs.thetaX, self.magnet.cs.thetaY, self.magnet.cs.thetaZ], degrees=True)
                rMagInv = rMag.inv()
                bMeasMag = rMagInv.apply(bMeasSys)
                print('bMeasMag:%s'%bMeasMag)
                mag2mgtrPosMag = self.magnet.updatePosition(bMeasMag[0],bMeasMag[1],bMeasMag[2],self.magnet.cs.x,self.magnet.cs.y,self.magnet.cs.z)
                mgtr.mag2mgtrPos = mag2mgtrPosMag
                #Below is vector pointing from magnetCS to magnetometerCS in
                mag2sysSys = rMag.apply(numpy.array([mag2mgtrPosMag[0],mag2mgtrPosMag[1],mag2mgtrPosMag[2]],numpy.double)) + [mgtr.cs.x,mgtr.cs.y,mgtr.cs.z]

                print('mag2sysSys:%s'%mag2sysSys)


            orientationError = 0
            for mgtr in self.magnetometers:
                mgtrActualPos = [mgtr.cs.x,mgtr.cs.y,mgtr.cs.z]
                for otherMgtr in self.magnetometers:
                    if otherMgtr != mgtr:
                        otherMgtrActualPos = [otherMgtr.cs.x, otherMgtr.cs.y, otherMgtr.cs.z]
                        actualDistance = [mgtrActualPos[0]-otherMgtrActualPos[0],mgtrActualPos[1]-otherMgtrActualPos[1],mgtrActualPos[2]-otherMgtrActualPos[2]]

                        measuredDistance = [-mgtr.mag2mgtrPos[0] + otherMgtr.mag2mgtrPos[0],
                                          -mgtr.mag2mgtrPos[1] + otherMgtr.mag2mgtrPos[1],
                                          -mgtr.mag2mgtrPos[2] + otherMgtr.mag2mgtrPos[2]]

                        errorDistance = [actualDistance[0]-measuredDistance[0],actualDistance[1]-measuredDistance[1],actualDistance[2]-measuredDistance[2]]

                        print('actualDistance: %s'%actualDistance)
                        print('measuredDistance: %s'%measuredDistance)
                        print('errorDistance: %s'%errorDistance)
                        orientationError+=math.sqrt(errorDistance[0]**2+errorDistance[1]**2+errorDistance[2]**2)

            print ('orientationError:%s'%orientationError)






            #print('Bmeas:%s'%Bmeas)

