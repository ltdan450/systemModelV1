import MHDCoordSys
import MHDRectPM
import MHDMagAssby
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
        self.magAssbys = []
        self.verbose = 1

    def addSystem(self,sys):
        self.systems.append(sys)

    def addMagnet(self,mag):
        self.magnets.append(mag)

    def addMagAssby(self, magAssby):
        self.magAssbys.append(magAssby)

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

    def getBFieldForA(self, a_sys, sys):
        # do da thang
        rotSys = R.from_euler(sys.cs.rot, [sys.cs.thetaX, sys.cs.thetaY, sys.cs.thetaZ], degrees=True)
        rInvSys = rotSys.inv()
        a_univ = rotSys.apply(a_sys)
        e_univ = numpy.array([sys.cs.x, sys.cs.y, sys.cs.z], numpy.double)
        c_univ = a_univ + e_univ
        if self.verbose:
            print('a_sys:%s \n a_univ: %s \n e_univ: %s \n c_univ: %s' % (a_sys, a_univ, e_univ, c_univ))

        bSum = numpy.array([0.0, 0.0, 0.0], numpy.double)
        for magAssby in self.magAssbys:
            d_univ = numpy.array([magAssby.cs.x, magAssby.cs.y, magAssby.cs.z], numpy.double)
            i_univ = c_univ - d_univ
            rotMagAssby =  R.from_euler(magAssby.cs.rot, [magAssby.cs.thetaX, magAssby.cs.thetaY, magAssby.cs.thetaZ], degrees=True)
            rInvMagAssby = rotMagAssby.inv()
            i_magAssby = rInvMagAssby.apply(i_univ)
            b_MagAssby = magAssby.getBForI(i_magAssby)
            b_univ = rotMagAssby.apply(b_MagAssby)
            b_sys = rInvSys.apply(b_univ)
            bSum += b_sys
            if self.verbose:
                print('d_univ:%s \n i_univ: %s \n i_magassby: %s \n b_magassby: %s' % (d_univ, i_univ, i_magAssby, b_MagAssby))


        return bSum





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
tests = [ 3]


if 1 in tests:
    # testing calculated magnetic fields vs measured ones obtain from:
        # https://www.andrews.edu/cas/physics/research/research_files/coilmag01.pdf
    cal_pts = {0.007:0.0341,0.016:0.0247,0.0258:0.0198,0.0359:0.0164,0.0459:0.0136,0.056:0.0111,0.0659:0.009,0.0759:0.0073,0.0858:0.006,0.096:0.0049,0.1059:0.004,0.1159:0.0033,0.1258:0.0028,0.1358:0.0023,0.1459:0.0019,0.1558:0.0016,0.1658:0.0015,0.1757:0.0013,0.1858:0.0011}

    xs = []
    bx_calcs = []
    bx_meass = []
    b_calcs = []
    xs_sorted = sorted(cal_pts.keys())

    for xcal in xs_sorted:
        u = MHDUniverse(.001, 130.0)

        xM = xcal
        yM = 0.0
        zM = 0.0
        magAssbyCS = MHDCoordSys.MHDCoordSys(xM, yM, zM, 0.0, 0.0, 0.0, u.cs)
        magAssby = MHDMagAssby.MHDMagAssby(magAssbyCS)
        u.addMagAssby(magAssby)

        magCS = MHDCoordSys.MHDCoordSys(0.0, 0.0, 0.0, 90.0, 000.0, 90.0, magAssbyCS)
        a = c.test_kingman_small['a']
        b = c.test_kingman_small['b']
        h = c.test_kingman_small['h']
        Br = c.test_kingman_small['Br']
        mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.0)
        mag.setJandKForBr(Br)
        magAssby.addMagnet(mag)

        #mgtrs = [[.0, .0, .0, 1]]  # , [.05, .05, .0, 0], [-.05, -.05, .0, 0]]
        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemB(csSys1, u)
        #sys1.magnet = mag
        u.addSystem(sys1)

        a_sys = numpy.array([0.0, 0.0, 0.0], numpy.double)

        bres = u.getBFieldForA(a_sys, sys1)
        #print('bres:%s'%bres)
        xs.append(xcal)
        bx_meass.append(cal_pts[xcal])
        b_calcs.append(bres)
        bx_calcs.append(bres[0])


    print('Comparing computed magnetic field to measured on p28 in coilmag01.pdf')
    for i in range(0,len(xs)):
        print('x:%s bxErrorPct:%s bxMeas:%s bxCalc:%s bCalc:%s'%(xs[i], round((bx_meass[i]-bx_calcs[i])/bx_meass[i]*100,2) ,bx_meass[i],bx_calcs[i],b_calcs[i]))
    #sys1.measInterval = 0.1

if 2 in tests:
    # testing calculated magnetic fields vs measured ones obtain from:
        # https://www.andrews.edu/cas/physics/research/research_files/coilmag01.pdf
    cal_pts = {0.081:0.1222,0.0867:0.0424,0.0968:0.0191,0.1069:0.0111,0.117:0.0073,0.1269:0.0051,0.1369:0.0036,0.1469:0.0028,0.1569:0.0023,0.167:0.0017,0.1769:0.0015,0.1868:0.0011,0.1969:0.0009,0.207:0.0008,0.2169:0.0008,0.227:0.0006,0.2369:0.0006,0.247:0.0006,0.2569:0.0006}

    ys = []
    bx_calcs = []
    bx_meass = []
    b_calcs = []
    ys_sorted = sorted(cal_pts.keys())

    for ycal in ys_sorted:
        u = MHDUniverse(.001, 130.0)

        xM = 0.0
        yM = ycal
        zM = 0.0
        magAssbyCS = MHDCoordSys.MHDCoordSys(xM, yM, zM, 0.0, 0.0, 0.0, u.cs)
        magAssby = MHDMagAssby.MHDMagAssby(magAssbyCS)
        u.addMagAssby(magAssby)

        magCS = MHDCoordSys.MHDCoordSys(0.0, 0.0, 0.0, 90.0, 000.0, 90.0, magAssbyCS)
        a = c.test_kingman_small['a']
        b = c.test_kingman_small['b']
        h = c.test_kingman_small['h']
        Br = c.test_kingman_small['Br']
        mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.0)
        mag.setJandKForBr(Br)
        magAssby.addMagnet(mag)

        #mgtrs = [[.0, .0, .0, 1]]  # , [.05, .05, .0, 0], [-.05, -.05, .0, 0]]
        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemB(csSys1, u)
        #sys1.magnet = mag
        u.addSystem(sys1)

        a_sys = numpy.array([0.0, 0.0, 0.0], numpy.double)

        bres = u.getBFieldForA(a_sys, sys1)
        #print('bres:%s'%bres)
        ys.append(ycal)
        bx_meass.append(-cal_pts[ycal])
        b_calcs.append(bres)
        bx_calcs.append(bres[0])


    print('Comparing computed magnetic field to measured on p29 in coilmag01.pdf')
    for i in range(0,len(ys)):
        print('y:%s bxErrorPct:%s bxMeas:%s bxCalc:%s bCalc:%s'%(ys[i], round((bx_meass[i]-bx_calcs[i])/bx_meass[i]*100,2) ,bx_meass[i],bx_calcs[i],b_calcs[i]))
    #sys1.measInterval = 0.1

if 3 in tests:
    # testing calculated magnetic fields vs measured ones obtain from:
        # https://www.andrews.edu/cas/physics/research/research_files/coilmag01.pdf
    cal_pts = {0.0549: 0.1277, 0.0608: 0.0403, 0.0709: 0.018, 0.081: 0.0102, 0.091: 0.0063, 0.1009: 0.0042,
                 0.111: 0.0032, 0.1211: 0.0023, 0.1311: 0.0019, 0.1411: 0.0014, 0.1511: 0.0012, 0.1611: 0.0009,
                 0.1711: 0.0009, 0.1811: 0.0007, 0.1912: 0.0005, 0.2012: 0.0005, 0.2112: 0.0004, 0.2212: 0.0005,
                 0.2312: 0.0004}

    zs = []
    bx_calcs = []
    bx_meass = []
    b_calcs = []
    zs_sorted = sorted(cal_pts.keys())

    for zcal in zs_sorted:
        u = MHDUniverse(.001, 130.0)

        xM = 0.0
        yM = 0.0
        zM = zcal
        magAssbyCS = MHDCoordSys.MHDCoordSys(xM, yM, zM, 0.0, 0.0, 0.0, u.cs)
        magAssby = MHDMagAssby.MHDMagAssby(magAssbyCS)
        u.addMagAssby(magAssby)

        magCS = MHDCoordSys.MHDCoordSys(0.0, 0.0, 0.0, 90.0, 000.0, 90.0, magAssbyCS)
        a = c.test_kingman_small['a']
        b = c.test_kingman_small['b']
        h = c.test_kingman_small['h']
        Br = c.test_kingman_small['Br']
        mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.0)
        mag.setJandKForBr(Br)
        magAssby.addMagnet(mag)

        #mgtrs = [[.0, .0, .0, 1]]  # , [.05, .05, .0, 0], [-.05, -.05, .0, 0]]
        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemB(csSys1, u)
        #sys1.magnet = mag
        u.addSystem(sys1)

        a_sys = numpy.array([0.0, 0.0, 0.0], numpy.double)

        bres = u.getBFieldForA(a_sys, sys1)
        #print('bres:%s'%bres)
        zs.append(zcal)
        bx_meass.append(-cal_pts[zcal])
        b_calcs.append(bres)
        bx_calcs.append(bres[0])


    print('Comparing computed magnetic field to measured on p30 in coilmag01.pdf')
    for i in range(0,len(zs)):
        print('z:%s bxErrorPct:%s bxMeas:%s bxCalc:%s bCalc:%s'%(zs[i], round((bx_meass[i]-bx_calcs[i])/bx_meass[i]*100,2) ,bx_meass[i],bx_calcs[i],b_calcs[i]))
    #sys1.measInterval = 0.1


if 4 in tests:
    # testing calculated magnetic fields vs measured ones obtain from:
        # https://www.andrews.edu/cas/physics/research/research_files/coilmag01.pdf
    cal_pts = {0.007:0.0341,0.016:0.0247,0.0258:0.0198,0.0359:0.0164,0.0459:0.0136,0.056:0.0111,0.0659:0.009,0.0759:0.0073,0.0858:0.006,0.096:0.0049,0.1059:0.004,0.1159:0.0033,0.1258:0.0028,0.1358:0.0023,0.1459:0.0019,0.1558:0.0016,0.1658:0.0015,0.1757:0.0013,0.1858:0.0011}

    xs = []
    bmagnitude_calcs = []
    bx_meass = []
    b_calcs = []
    xs_sorted = sorted(cal_pts.keys())

    for xcal in xs_sorted:
        u = MHDUniverse(.001, 130.0)

        a = xcal/math.sqrt(2.0)

        xM = a
        yM = a
        zM = 0.0
        magAssbyCS = MHDCoordSys.MHDCoordSys(xM, yM, zM, 0.0, 0.0, 45.0, u.cs)
        magAssby = MHDMagAssby.MHDMagAssby(magAssbyCS)
        u.addMagAssby(magAssby)

        magCS = MHDCoordSys.MHDCoordSys(0.0, 0.0, 0.0, 90.0, 000.0, 90.0, magAssbyCS)
        a = c.test_kingman_small['a']
        b = c.test_kingman_small['b']
        h = c.test_kingman_small['h']
        Br = c.test_kingman_small['Br']
        mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.0)
        mag.setJandKForBr(Br)
        magAssby.addMagnet(mag)

        #mgtrs = [[.0, .0, .0, 1]]  # , [.05, .05, .0, 0], [-.05, -.05, .0, 0]]
        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemB(csSys1, u)
        #sys1.magnet = mag
        u.addSystem(sys1)

        a_sys = numpy.array([0.0, 0.0, 0.0], numpy.double)

        bres = u.getBFieldForA(a_sys, sys1)
        #print('bres:%s'%bres)
        xs.append(xcal)
        bx_meass.append(cal_pts[xcal])
        b_calcs.append(bres)
        bmagnitude_calcs.append(math.sqrt(bres[0]*bres[0]+bres[1]*bres[1]+bres[2]*bres[2]))


    print('Computing magnetic field for magnet rotated 45deg and moved 1/sqrt(2)*a in x and y'
          '\n where a is the length along the z axis in shown on p28 in coilmag01.pdf'
          '\n the magnitude of the magnetic field is compared with the measured z axis values on that page')
    for i in range(0,len(xs)):
        print('x:%s bxErrorPct:%s bxMeas:%s bmagCalc:%s bCalc:%s'%(xs[i], round((bx_meass[i]-bmagnitude_calcs[i])/bx_meass[i]*100,2) ,bx_meass[i],bmagnitude_calcs[i],b_calcs[i]))
    #sys1.measInterval = 0.1

if 5 in tests:
    # testing calculated magnetic fields vs measured ones obtain from:
        # https://www.andrews.edu/cas/physics/research/research_files/coilmag01.pdf
    cal_pts = {0.007:0.0341,0.016:0.0247,0.0258:0.0198,0.0359:0.0164,0.0459:0.0136,0.056:0.0111,0.0659:0.009,0.0759:0.0073,0.0858:0.006,0.096:0.0049,0.1059:0.004,0.1159:0.0033,0.1258:0.0028,0.1358:0.0023,0.1459:0.0019,0.1558:0.0016,0.1658:0.0015,0.1757:0.0013,0.1858:0.0011}

    xs = []
    bx_calcs = []
    bx_meass = []
    b_calcs = []
    xs_sorted = sorted(cal_pts.keys())

    for xcal in xs_sorted:
        # setup universe
        u = MHDUniverse(.001, 130.0)

        # setup mag assby
        xM = xcal
        yM = 0.0
        zM = 0.0
        magAssbyCS = MHDCoordSys.MHDCoordSys(xM, yM, zM, 0.0, 0.0, 0.0, u.cs)
        magAssby = MHDMagAssby.MHDMagAssby(magAssbyCS)
        u.addMagAssby(magAssby)
        magCS = MHDCoordSys.MHDCoordSys(0.0, 0.0, 0.0, 90.0, 000.0, 90.0, magAssbyCS)
        a = c.test_kingman_small['a']
        b = c.test_kingman_small['b']
        h = c.test_kingman_small['h']
        Br = c.test_kingman_small['Br']
        mag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.0)
        mag.setJandKForBr(Br)
        magAssby.addMagnet(mag)

        #setup system
        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .0, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemB(csSys1, u)
        u.addSystem(sys1)
        mgtr = MHDMagnetometer.MHDMagnetometer(MHDCoordSys.MHDCoordSys(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, sys1.cs))
        sys1.addMagnetometer(mgtr)
        #setup virtual magnet for system
        vMagAssbyCS = MHDCoordSys.MHDCoordSys(xM, yM, zM, 0.0, 0.0, 0.0, u.cs)
        vMagAssby = MHDMagAssby.MHDMagAssby(vMagAssbyCS)
        vMagCS = MHDCoordSys.MHDCoordSys(0.0, 0.0, 0.0, 90.0, 000.0, 90.0, vMagAssbyCS)
        vMag = MHDRectPM.MHDRectPM(vMagCS, a, b, h, 1.0)
        vMag.setJandKForBr(Br)
        vMagAssby.addMagnet(vMag)
        sys1.addMagAssby(vMagAssby)

        sys1.getPos()



        if 1 == 0:
            a_sys = numpy.array([0.0, 0.0, 0.0], numpy.double)
            bres = u.getBFieldForA(a_sys, sys1)
            #print('bres:%s'%bres)
            xs.append(xcal)
            bx_meass.append(cal_pts[xcal])
            b_calcs.append(bres)
            bx_calcs.append(bres[0])


    print('Comparing computed magnetic field to measured on p28 in coilmag01.pdf')
    for i in range(0,len(xs)):
        print('x:%s bxErrorPct:%s bxMeas:%s bxCalc:%s bCalc:%s'%(xs[i], round((bx_meass[i]-bx_calcs[i])/bx_meass[i]*100,2) ,bx_meass[i],bx_calcs[i],b_calcs[i]))
    #sys1.measInterval = 0.1


if 33 in tests:
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

if -44 in tests:
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
        a = 1.5     * .0254  # 0.048
        b = 1.0     * .0254  # 0.022
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
        paramset = [{'noiseTesla': 1.0e-22, 'sensitivityTesla': 1.0e-22},{'noiseTesla': 3.0e-7, 'sensitivityTesla': 1.0e-8},{'noiseTesla': 3.8e-8, 'sensitivityTesla': 1.3e-8}]
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