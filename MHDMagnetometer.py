import MHDCoordSys
import MHDRectPM
import numpy as np
from scipy.spatial.transform import Rotation as R
import MHDTools
from scipy import signal
from matplotlib import pyplot as plt



class MHDMagnetometer:

    def __init__(self,cs=0,system=0,paramsin={}):

        try:
            self.noiseSigma = paramsin['noiseTesla']
        except:
            self.noiseSigma = 0.0
        try:
            self.sensitivity = paramsin['sensitivityTesla']
        except:
            self.sensitivity = 1.0e-12
        try:
            self.maxB = paramsin['maxBTesla']
        except:
            self.maxB = 1.0e12
        try:
            self.cost = paramsin['cost']
        except:
            self.cost = 1.0e6


        self.cs = cs
        self.altCS = 0
        self.system = system
        if system!=0:
            self.universe = system.universe
        self.mag2mgtrPos = [0.0, 0.0, 0.0]
        self.noise = []
        self.configure()

        self.filtStarted = 0


    def configure(self,bits=16,sensitivityLSBPerG=12000):
        #sensitivity of QMC5883L is 12000 LSB/G
        doNothing = 1
        #bitsPerG = 12000
        #bitsPerT = bitsPerG * 10000
        #tPerBit = 1.0/bitsPerT

        #self.sensitivity = tPerBit
    def configFilt(self):
        if 1 == 0:
            fs = 200
            fc = 5

            filt = signal.iirdesign(wp=0.01, ws=0.05, gstop=60, gpass=1, output='ba',ftype='ellip')
            #print('filt:%s'%filt)
            self.b = filt[0]
            self.a = filt[1]
            #self.b = signal.firwin(90,.02)
            #self.a = 1

            self.zX = signal.lfilter_zi(self.b, self.a)
            self.zY = signal.lfilter_zi(self.b, self.a)
            self.zZ = signal.lfilter_zi(self.b, self.a)

        if 1 == 1:
            fs = 100
            N = 10
            #scipy.signal.iirfilter(N, Wn, rp=None, rs=None, btype='band', analog=False, ftype='butter', output='ba', fs=None)[source]
            self.soss = signal.iirfilter(N, .05, btype = 'lowpass', analog = False, ftype = 'butter', output='sos')
            self.zX = self.sosfilt_zi(self.soss)
            self.zY = self.sosfilt_zi(self.soss)
            self.zZ = self.sosfilt_zi(self.soss)
            if 1 == 0:
                w, h = signal.sosfreqz(self.soss)
                #print('w:%s h:%s'%(w,h))
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax.semilogx(w / (2 * np.pi), 20 * np.log10(np.maximum(abs(h), 1e-5)))
                ax.set_title('Chebyshev Type II bandpass frequency response')
                ax.set_xlabel('Frequency [Hz]')
                ax.set_ylabel('Amplitude [dB]')
                ax.axis((.001, 1, -100, 10))
                ax.grid(which='both', axis='both')
                #plt.show()

                #self.a = a
                #self.b = b

                #self.zX = signal.lfilter_zi(self.b, self.a)
                #self.zY = signal.lfilter_zi(self.b, self.a)
                #self.zZ = signal.lfilter_zi(self.b, self.a)
    def filterData(self,data):

        filtType = 'not sos'

        if not self.filtStarted:
            for i in range (0,100,1):
                if filtType == 'sos':
                    filtX, self.zX = signal.sosfilt(self.soss,[data[0]],zi=self.zX)#signal.lfilter(self.b, self.a, [data[0]], zi=self.zX)
                    filtY, self.zY = signal.sosfilt(self.soss,[data[1]],zi=self.zY)#signal.lfilter(self.b, self.a, [data[1]], zi=self.zY)
                    filtZ, self.zZ = signal.sosfilt(self.soss,[data[2]],zi=self.zZ)#signal.lfilter(self.b, self.a, [data[2]], zi=self.zZ)
                else:
                    filtX, self.zX = signal.lfilter(self.b, self.a, [data[0]], zi=self.zX)
                    filtY, self.zY = signal.lfilter(self.b, self.a, [data[1]], zi=self.zY)
                    filtZ, self.zZ = signal.lfilter(self.b, self.a, [data[2]], zi=self.zZ)
            self.filtStarted = 1

        if filtType == 'sos':
            filtX, self.zX = signal.sosfilt(self.soss, [data[0]], zi=self.zX)
            filtY, self.zY = signal.sosfilt(self.soss, [data[1]], zi=self.zY)
            filtZ, self.zZ = signal.sosfilt(self.soss, [data[2]], zi=self.zZ)
        else:
            filtX, self.zX = signal.lfilter(self.b, self.a, [data[0]], zi=self.zX)
            filtY, self.zY = signal.lfilter(self.b, self.a, [data[1]], zi=self.zY)
            filtZ, self.zZ = signal.lfilter(self.b, self.a, [data[2]], zi=self.zZ)


        out = [filtX[0], filtY[0], filtZ[0]]
        #print('filt:%s data:%s'%(out,data))


        return out
    def sosfilt_zi(self,sos):
        """Compute an initial state `zi` for the sosfilt function"""
        from scipy.signal import lfilter_zi
        sos = np.asarray(sos)
        if sos.ndim != 2 or sos.shape[1] != 6:
            raise ValueError('sos must be shape (n_sections, 6)')

        n_sections = sos.shape[0]
        zi = np.empty((n_sections, 2))
        scale = 1.0
        for section in range(n_sections):
            b = sos[section, :3]
            a = sos[section, 3:]
            zi[section] = scale * lfilter_zi(b, a)
            # If H(z) = B(z)/A(z) is this section's transfer function, then
            # b.sum()/a.sum() is H(1), the gain at omega=0.  That's the steady
            # state value of this section's step response.
            scale *= b.sum() / a.sum()

        return zi

    ######### normal operation methods#########
    def measure(self):
        #print('measure start')
        bFieldIn = self.getBMeasurementFromUniverse()
        #print(bFieldIn)
        #print('sens:%s  bFieldIn2:%s'%(self.sensitivity,bFieldIn[2]))
        #print('intsZ:%s'%(round(bFieldIn[2]/self.sensitivity)*self.sensitivity))

        if 1 == 1:
            bFieldIn[0] = round(bFieldIn[0] / self.sensitivity) * self.sensitivity + self.getNoise()
            bFieldIn[1] = round(bFieldIn[1] / self.sensitivity) * self.sensitivity + self.getNoise()
            bFieldIn[2] = round(bFieldIn[2] / self.sensitivity) * self.sensitivity + self.getNoise()

        if self.altCS and 1 == 0:
            rot = R.from_euler('xyz',[self.altCS.thetaX,self.altCS.thetaY,self.altCS.thetaZ],degrees=True)
            #bFieldIn = rot.inv().apply(bFieldIn)

        if 1 == 1:
            bFieldIn = self.filterData(bFieldIn)

        return bFieldIn

    #########omnipotent methods#########
    def getBMeasurementFromUniverse(self):
        bField = self.system.universe.getBFieldForMgtr(self)

        #print('bField:%s'%bField)
        return bField

    def getNoise(self):
        if len(self.noise)==0:
            self.noise = np.random.normal(0,self.noiseSigma+1e-11,1000)

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
    mag = MHDMagnetometer(0,0,0)
    mag.testNoise(1.0/fs,t)