import MHDCoordSys
import numpy
from scipy.spatial.transform import Rotation as R
import MHDTools



class MHDMagnetometer:

    def __init__(self,cs=0,system=0):
        self.cs = cs
        self.altCS = 0
        self.system = system
        if system!=0:
            self.universe = system.universe
        self.mag2mgtrPos = [0.0, 0.0, 0.0]
        self.noise = []

    def configure(self,bits=16,sensitivity=12000):
        self.sensitivity = sensitivity

    ######### normal operation methods#########
    def measure(self):
        #print('measure start')
        bFieldIn = self.getBMeasurementFromUniverse()

        bFieldIn[0] = round(bFieldIn[0] / self.sensitivity) * self.sensitivity
        bFieldIn[1] = round(bFieldIn[1] / self.sensitivity) * self.sensitivity
        bFieldIn[2] = round(bFieldIn[2] / self.sensitivity) * self.sensitivity

        if self.altCS:
            rot = R.from_euler('xyz',[self.altCS.thetaX,self.altCS.thetaY,self.altCS.thetaZ],degrees=True)
            bFieldIn = rot.inv().apply(bFieldIn)

        return bFieldIn

    #########omnipotent methods#########
    def getBMeasurementFromUniverse(self):
        bField = self.system.universe.getBFieldForMgtr(self)
        #print('bField:%s'%bField)
        return bField

    def getNoise(self,stdev=0.0025):
        if len(self.noise)==0:
            self.noise = numpy.random.normal(0,stdev,1000)

        last, self.noise = self.noise[-1], self.noise[:-1]
        return last


    def testNoise(self,dt=1/80,tMax=1):
        ts = []
        bs = []
        for i in range(0,int(tMax/dt),1):
            ts.append(dt * i)
            bs.append(self.getNoise(0.0025))

        MHDTools.plotPSD(ts,bs)

if 1 == 0:
    fs = 80.0
    nSamples = 10000
    t = nSamples/fs
    mag = MHDMagnetometer(0,0)
    mag.testNoise(1.0/fs,t)