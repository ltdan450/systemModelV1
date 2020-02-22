import MHDCoordSys
import MHDRectPM
import MHDSystem
import MHDMagnetometer
import numpy
import math
import scipy.signal as signal
from scipy.spatial.transform import Rotation as R
import globalConstants as c

def log(toLog):
    verbose = False
    if verbose:
        print(toLog)

def sosfilt_zi(sos):
    """Compute an initial state `zi` for the sosfilt function"""
    from scipy.signal import lfilter_zi
    sos = numpy.asarray(sos)
    if sos.ndim != 2 or sos.shape[1] != 6:
        raise ValueError('sos must be shape (n_sections, 6)')

    n_sections = sos.shape[0]
    zi = numpy.empty((n_sections, 2))
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

class MHDUniverse:

    def __init__(self,dt,tMax):
        self.dt = dt
        self.tMax = tMax
        self.cs = MHDCoordSys.MHDCoordSys()
        self.systems = []
        self.magnets = []
        self.mainMag = 0
        self.rMagZero = [0.0,0.0,0.0]


    def addSystem(self,sys):
        self.systems.append(sys)

    def addMagnet(self,mag):
        self.magnets.append(mag)


    def getBFieldForMgtr(self,mgtr):
        #print('getBFieldForMgtr')
        #print('\n')
        bMgtrSum = numpy.array([0.0, 0.0, 0.0], numpy.double)
        for mag in self.magnets:
            #print('eval mag:%s'%mag)
            sys = mgtr.system
            #aMgtr = numpy.array([mgtr.cs.x,mgtr.cs.y,mgtr.cs.z],numpy.double)
            #rMgtr = R.from_euler(mgtr.cs.rot, [mgtr.cs.thetaX, mgtr.cs.thetaY, mgtr.cs.thetaZ], degrees=True)
            aSys = numpy.array([mgtr.cs.x,mgtr.cs.y,mgtr.cs.z],numpy.double)
            rotSys = R.from_euler(sys.cs.rot, [sys.cs.thetaX, sys.cs.thetaY, sys.cs.thetaZ], degrees=True)
            aUniv = rotSys.apply(aSys)
            eSys = numpy.array([sys.cs.x,sys.cs.y,sys.cs.z],numpy.double)
            eUniv = rotSys.apply(eSys)
            cUniv = eUniv + aUniv
            #log('rMag = %s'%)

            #dMag = numpy.array([mag.cs.x, mag.cs.y, mag.cs.z], numpy.double)
            rotMag = R.from_euler(mag.cs.rot, [mag.cs.thetaX, mag.cs.thetaY, mag.cs.thetaZ], degrees=True)

            dUniv = numpy.array([mag.cs.x, mag.cs.y, mag.cs.z], numpy.double)
            fUniv = cUniv - dUniv

            rInvMag = rotMag.inv()
            fMag = rInvMag.apply(fUniv)

            log('f_univ:%s \nf_mag:%s'%(fUniv,fMag))

            #need to get b field now from fMag
            bMag = mag.getB(fMag)
            bUniv = rotMag.apply(bMag)
            rSysInv = rotSys.inv()
            bSys = rSysInv.apply(bUniv)

            rMgtr = R.from_euler(mgtr.cs.rot, [mgtr.cs.thetaX, mgtr.cs.thetaY, mgtr.cs.thetaZ], degrees=True)
            rMgtrInv = rMgtr.inv()
            bMgtr = rMgtrInv.apply(bSys)
            log('b_mag:%s \nb_mgtr:%s'%(bMag,bMgtr))
            bMgtrSum+= bMgtr

        return bMgtrSum

    def run(self):
        t = 0
        rateErrors = []
        ampErrors = []

        while t<=self.tMax:
            #print('t:%s'%(t))
            #update positions

            mag = .002
            f = .5
            tGrace = 20

            self.mainMag.cs.z = self.rMagZero[2] + 1* mag * math.sin(f*t*2*math.pi) + 0* .2* mag * math.sin(20*t*2*math.pi)

            #self.mainMag.cs.thetaY = 0 + 1 * math.sin(f * t * 2 * math.pi)
            #print('ty:%s'%(self.mainMag.cs.thetaY))
            #update system
            for system in self.systems:
                res = system.step(t)
                if res:
                    rate,amp = res
                    if t>tGrace:
                        rateErrors.append(abs(f-rate))
                        ampErrors.append(abs(amp - (2*mag)))


            t+=self.dt
        for system in self.systems:
            system.testEnd()
        avgRateError = sum(rateErrors)/len(rateErrors)
        avgAmpError = sum(ampErrors) / len(ampErrors)
        print('i:%s avgRateError:%s avgAmpError:%s'%(i,avgRateError,avgAmpError))


tests = c.tests
tests = [-6]

if 3 in tests:
    for i1 in range(0,1,1):
        #setup universe
        u = MHDUniverse(.001,10)

        #setup magnets
        a = 1.2*.0254#0.048
        b = 0.6*.0254#0.022
        h = 0.1875*.0254#0.011
        magCS = MHDCoordSys.MHDCoordSys(0.0, -0.00, -0.3, 0, -0, 90, u.cs)
        mag = MHDRectPM.MHDRectPM(magCS,a, b, h, 1.0)
        mag.setJandKForBr(1.31)
        u.addMagnet(mag)
        u.mainMag = mag

        #setup system
        mgtrs = [[0,0,0,0],[0.05,0,0,0],[-0.03,.01,0,0]]

        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemA(csSys1, u)
        sys1.magnet=mag
        u.addSystem(sys1)

        #setup magnetometers
        paramset = [{'noiseTesla': 1.0e-22, 'sensitivityTesla': 1.0e-22}]
        for mgtr in mgtrs:
            mgtrCS = MHDCoordSys.MHDCoordSys(mgtr[0],mgtr[1],mgtr[2],0,0,0,sys1.cs)
            sys1.addMagnetometer(MHDMagnetometer.MHDMagnetometer(mgtrCS,sys1,paramset[mgtr[3]]))

        if 1 == 1:
            #test magnet
            if 1==0:
                for i in range(-10,10,1):
                    testMagCS = MHDCoordSys.MHDCoordSys(0, 0, 0.0, 00.0, 0.0, 0.0, u.cs)
                    testMag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.25e11)
                    dx = .025
                    pos = [0*i*dx,0*i*dx,0+i*dx]
                    print('testMagB:%s pos:%s'%(testMag.getB(pos),pos))

            dx = .2
            var = -1 + dx * i1
            #print('var:%s'%var)
            vMagCS = MHDCoordSys.MHDCoordSys(.0, .0, -0.270, 10.0, 10.0, 0.0, sys1.cs)
            vMag = MHDRectPM.MHDRectPM(vMagCS, a, b, h, 1.25e11*1)
            vMag.verbose = 1
            sys1.addMagnet(vMag)

            u.run()

if -4 in tests:
    u = MHDUniverse(.001, 30.0)

    # setup magnets
    a = c.a #1.2 * .0254  # 0.048
    b = c.b #0.6 * .0254  # 0.022
    h = c.h #0.1875 * .0254  # 0.011
    magCS = MHDCoordSys.MHDCoordSys(c.xRm, c.yRm, c.zRm, c.tXRm, c.tYRm, c.tZRm, u.cs)
    mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.0)
    mag.setJandKForBr(c.Br)
    u.addMagnet(mag)
    u.mainMag = mag

    # setup system
    mgtrs = [[c.x1, c.y1, c.z1, 0], [c.x2, c.y2, c.z2, 0], [.04, c.y2/2, .06, 0]]
    csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
    sys1 = MHDSystem.MHDSystemA(csSys1, u)
    sys1.magnet = mag
    u.addSystem(sys1)

    sys1.measInterval = 0.05

    #

    #system high pass filter
    if 1 == 1:
        bpmThresh = 10
        fThresh = bpmThresh / 60.0
        fHigh = 2. / 60.
        fLow = 1. / 60.
        ord = 1
        fNyq = 0.5 * 1.0 / sys1.measInterval
        print('fH:%s fL:%s nyq:%s' % (fThresh, fLow, fNyq))
        b, a = signal.butter(ord, [fThresh / fNyq], btype='highpass')
        sys1.posFilt = [b, a]
        sys1.filtType='FIR'

    #system low pass filter
    if 1 == 1:
        bpmThresh = 200
        fThresh = bpmThresh / 60.0
        ord = 4
        fNyq = 0.5 * 1.0 / sys1.measInterval
        #print('fH:%s fL:%s nyq:%s' % (fThresh, fLow, fNyq))
        b, a = signal.butter(ord, [fThresh / fNyq], btype='lowpass')
        lowpassFilt = [b, a]
        #sys1.filtType='FIR'


    # setup magnetometers
    paramset = [{'noiseTesla': 1.0e-22, 'sensitivityTesla': 1.0e-22}]
    for mgtr in mgtrs:
        mgtrCS = MHDCoordSys.MHDCoordSys(mgtr[0], mgtr[1], mgtr[2], 0, 0, 0, sys1.cs)
        newMgtr = MHDMagnetometer.MHDMagnetometer(mgtrCS, sys1, paramset[mgtr[3]])
        sys1.addMagnetometer(newMgtr)
        newMgtr.a = a
        newMgtr.b = b
        if 1 == 1:
            newMgtr.zX = signal.lfilter_zi(newMgtr.b, newMgtr.a)
            newMgtr.zY = signal.lfilter_zi(newMgtr.b, newMgtr.a)
            newMgtr.zZ = signal.lfilter_zi(newMgtr.b, newMgtr.a)


    u.run()

# sensor array under mattress
if -5 in tests:

    for i in range (0,24,4):
        u = MHDUniverse(.001, 130.0)

        # setup magnets
        a = 1.2     * .0254  # 0.048
        b = 0.6     * .0254  # 0.022
        h = 0.1875  * .0254  # 0.011

        xIn = 9.4 * 0 + 1.*i *0.
        yIn = 13.6 * 0.
        zIn = 1. * i
        xM = xIn*2.54/100.
        yM = yIn*2.54/100.
        zM = zIn*2.54/100.
        print('xM:%s yM:%s zM:%s'%(xM,yM,zM))

        magCS = MHDCoordSys.MHDCoordSys(xM, yM, zM, .0, .0, .0, u.cs)
        mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.0)
        u.rMagZero = [magCS.x,magCS.y,magCS.z]
        mag.setJandKForBr(1.3)
        u.addMagnet(mag)
        u.mainMag = mag

        # setup system
        mgtrs = [[.0, .0, .0, 1]]#, [.05, .05, .0, 0], [-.05, -.05, .0, 0]]
        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemA(csSys1, u)
        sys1.magnet = mag
        u.addSystem(sys1)

        sys1.measInterval = 0.1

        #system high pass filter
        if 1 == 1:
            bpmLowThresh = 20.
            bpmHighThresh = 65.
            fLow= bpmLowThresh / 60.0
            fHigh = bpmHighThresh / 60.
            ord = 3
            fNyq = 0.5 * 1.0 / sys1.measInterval
            #print('fH:%s fL:%s nyq:%s' % (bpmHighThresh, bpmLowThresh, fNyq))
            #b, a = signal.butter(ord, [fLow / fNyq], btype='highpass')
            #b, a = signal.butter(ord, [fHigh/fNyq,fLow / fNyq], btype='bandpass')
            b, a = signal.butter(ord, [ fLow / fNyq, fHigh / fNyq], btype='bandpass')
            sys1.posFilt = [b, a]
            sys1.filtType='FIR'

        #system low pass filter
        if 1 == 1:
            bpmThresh = 65.
            fThresh = bpmThresh / 60.0
            ord = 6
            fNyq = 0.5 * 1.0 / sys1.measInterval
            #print('fH:%s fL:%s nyq:%s' % (fThresh, fLow, fNyq))
            b, a = signal.butter(ord, [fThresh / fNyq], btype='lowpass')
            lowpassFilt = [b, a]
            #sys1.filtType='FIR'


        # setup magnetometers
        paramset = [{'noiseTesla': 1.0e-22, 'sensitivityTesla': 1.0e-22},{'noiseTesla': 3.0e-7, 'sensitivityTesla': 1.0e-8},{'noiseTesla': 1.5e-8, 'sensitivityTesla': 1.3e-8}]
        for mgtr in mgtrs:
            mgtrCS = MHDCoordSys.MHDCoordSys(mgtr[0], mgtr[1], mgtr[2], 0, 0, 0, sys1.cs)
            newMgtr = MHDMagnetometer.MHDMagnetometer(mgtrCS, sys1, paramset[mgtr[3]])
            sys1.addMagnetometer(newMgtr)
            newMgtr.a = a
            newMgtr.b = b
            if 1 == 1:
                newMgtr.zX = signal.lfilter_zi(newMgtr.b, newMgtr.a)
                newMgtr.zY = signal.lfilter_zi(newMgtr.b, newMgtr.a)
                newMgtr.zZ = signal.lfilter_zi(newMgtr.b, newMgtr.a)


        u.run()


if -6 in tests:

    for i in range (20,24, 4):
        u = MHDUniverse(.001, 50.0)

        # setup magnets
        a = 1.2     * .0254  # 0.048
        b = 0.6     * .0254  # 0.022
        h = 0.1875  * .0254  # 0.011

        xIn = 9.4 * 0 + 1.*i
        yIn = 11 * 1.
        zIn = 3
        xM = xIn*2.54/100.
        yM = yIn*2.54/100.
        zM = zIn*2.54/100.

        print('xM:%s yM:%s zM:%s'%(xM,yM,zM))

        magCS = MHDCoordSys.MHDCoordSys(xM, yM, zM, .0, .0, .0, u.cs)
        mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.0)
        u.rMagZero = [magCS.x,magCS.y,magCS.z]
        mag.setJandKForBr(1.3)
        u.addMagnet(mag)
        u.mainMag = mag

        # setup system
        mgtrs = [[.0, .0, .0, 2]]#, [.05, .05, .0, 0], [-.05, -.05, .0, 0]]
        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemA(csSys1, u)
        sys1.magnet = mag
        u.addSystem(sys1)

        sys1.measInterval = 0.1

        #system high pass filter
        if 1 == 1:
            bpmLowThresh = 20.
            bpmHighThresh = 65.
            fLow= bpmLowThresh / 60.0
            fHigh = bpmHighThresh / 60.
            ord = 3
            fNyq = 0.5 * 1.0 / sys1.measInterval
            #print('fH:%s fL:%s nyq:%s' % (bpmHighThresh, bpmLowThresh, fNyq))
            #b, a = signal.butter(ord, [fLow / fNyq], btype='highpass')
            #b, a = signal.butter(ord, [fHigh/fNyq,fLow / fNyq], btype='bandpass')
            b, a = signal.butter(ord, [ fLow / fNyq, fHigh / fNyq], btype='bandpass')
            sys1.posFilt = [b, a]
            sys1.filtType='FIR'

        #system low pass filter
        if 1 == 1:
            bpmThresh = bpmHighThresh
            fThresh = bpmThresh / 60.0
            ord = 6
            fNyq = 0.5 * 1.0 / sys1.measInterval
            #print('fH:%s fL:%s nyq:%s' % (fThresh, fLow, fNyq))
            b, a = signal.butter(ord, [fThresh / fNyq], btype='lowpass')
            lowpassFilt = [b, a]
            #sys1.filtType='FIR'


        # setup magnetometers
        paramset = [{'noiseTesla': 1.0e-22, 'sensitivityTesla': 1.0e-22},{'noiseTesla': 3.0e-7, 'sensitivityTesla': 1.0e-8},{'noiseTesla': 1.5e-8, 'sensitivityTesla': 1.3e-8}]
        for mgtr in mgtrs:
            mgtrCS = MHDCoordSys.MHDCoordSys(mgtr[0], mgtr[1], mgtr[2], 0, 0, 0, sys1.cs)
            newMgtr = MHDMagnetometer.MHDMagnetometer(mgtrCS, sys1, paramset[mgtr[3]])
            sys1.addMagnetometer(newMgtr)
            newMgtr.a = a
            newMgtr.b = b
            if 1 == 1:
                newMgtr.zX = signal.lfilter_zi(newMgtr.b, newMgtr.a)
                newMgtr.zY = signal.lfilter_zi(newMgtr.b, newMgtr.a)
                newMgtr.zZ = signal.lfilter_zi(newMgtr.b, newMgtr.a)


        u.run()



#just magnets and magnetometers
if 1 == 0:

    # setup universe coordinate system
    uCS = MHDCoordSys.MHDCoordSys()

    # setup magnet
    a = .047
    b = .0215
    h = 0.0097
    magCS = MHDCoordSys.MHDCoordSys(0.0, 0.0, 0.0, 0, 0, 0.0, u.cs)
    mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.25e11 * 1)



    mgtr1CS = MHDCoordSys.MHDCoordSys(.0, .0, 0.10, 0.0, 0.0, 0.0, u.cs)
    mgtr1 = MHDMagnetometer.MHDMagnetometer(mgtr1CS, 0)



    rMgtr1_mag1 = rotMgtr1_univ.inv().apply(rMgtr1_univ)

    mgtr2CS = MHDCoordSys.MHDCoordSys(.020, .0, 0.10, 0.0, 0.0, 0.0, u.cs)
    #mgtr2 = MHDMagnetometer.MHDMagnetometer(mgtr2CS, 0)











    u.addMagnet(mag)

    # setup magnetometers
    mgtr1CS = MHDCoordSys.MHDCoordSys(.0, .0, 0.10, 0.0, 0.0, 0.0, sys1.cs)
    mgtr1 = MHDMagnetometer.MHDMagnetometer(mgtr1CS, sys1)
    sys1.magnetometers.append(mgtr1)

    mgtr2CS = MHDCoordSys.MHDCoordSys(0.0200, .0, -.0, 0.0, 0.0, 0.0, sys1.cs)
    mgtr2 = MHDMagnetometer.MHDMagnetometer(mgtr2CS, sys1)
    # sys1.magnetometers.append(mgtr2)

    vMagCS = MHDCoordSys.MHDCoordSys(.0, .0, 0.270, 0.0, 0.0, 0.0, sys1.cs)
    vMag = MHDRectPM.MHDRectPM(vMagCS, a, b, h, 1.25e11 * 1)
    vMag.verbose = 1
    sys1.addMagnet(vMag)