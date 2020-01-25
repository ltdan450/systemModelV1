import MHDCoordSys
import numpy
from scipy.spatial.transform import Rotation as R
import MHDTools
from scipy import signal


class MHDMagnetometer:

    def __init__(self,cs=0,system=0):
        self.cs = cs
        self.altCS = 0
        self.system = system
        if system!=0:
            self.universe = system.universe
        self.mag2mgtrPos = [0.0, 0.0, 0.0]
        self.noise = []
        self.configure()
        self.configFilt()

        self.filtStarted = 0

    def configure(self,bits=16,sensitivityLSBPerG=12000):
        #sensitivity of QMC5883L is 12000 LSB/G
        bitsPerG = 12000
        bitsPerT = bitsPerG * 10000
        tPerBit = 1.0/bitsPerT

        self.sensitivity = tPerBit

    def configFilt(self):

        fs = 200
        fc = 5

        self.b = signal.firwin(90,.01)

        self.zX = signal.lfilter_zi(self.b, 1)
        self.zY = signal.lfilter_zi(self.b, 1)
        self.zZ = signal.lfilter_zi(self.b, 1)



        # for i, x in enumerate(data):
        #    result[i], z = signal.lfilter(b, 1, [x], zi=z)

    def filterData(self,data):

        if not self.filtStarted:
            for i in range (0,1000,1):
                filtX, self.zX = signal.lfilter(self.b, 1.0, [data[0]], zi=self.zX)
                filtY, self.zY = signal.lfilter(self.b, 1.0, [data[1]], zi=self.zY)
                filtZ, self.zZ = signal.lfilter(self.b, 1.0, [data[2]], zi=self.zZ)
        self.filtStarted = 1


        filtX, self.zX = signal.lfilter(self.b, 1.0, [data[0]], zi=self.zX)
        filtY, self.zY = signal.lfilter(self.b, 1.0, [data[1]], zi=self.zY)
        filtZ, self.zZ = signal.lfilter(self.b, 1.0, [data[2]], zi=self.zZ)

        out = [filtX[0], filtY[0], filtZ[0]]
        print('filt:%s'%(out))


        return out


    ######### normal operation methods#########
    def measure(self):
        #print('measure start')
        bFieldIn = self.getBMeasurementFromUniverse()
        #print(bFieldIn)
        #print('sens:%s  bFieldIn2:%s'%(self.sensitivity,bFieldIn[2]))
        #print('intsZ:%s'%(round(bFieldIn[2]/self.sensitivity)*self.sensitivity))

        bFieldIn[0] = round(bFieldIn[0] / self.sensitivity) * self.sensitivity + self.getNoise()
        bFieldIn[1] = round(bFieldIn[1] / self.sensitivity) * self.sensitivity + self.getNoise()
        bFieldIn[2] = round(bFieldIn[2] / self.sensitivity) * self.sensitivity + self.getNoise()

        if self.altCS:
            rot = R.from_euler('xyz',[self.altCS.thetaX,self.altCS.thetaY,self.altCS.thetaZ],degrees=True)
            bFieldIn = rot.inv().apply(bFieldIn)

        bFieldIn = self.filterData(bFieldIn)

        return bFieldIn

    #########omnipotent methods#########
    def getBMeasurementFromUniverse(self):
        bField = self.system.universe.getBFieldForMgtr(self)

        #print('bField:%s'%bField)
        return bField

    def getNoise(self,stdev=0.0025/10000.0):
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