import MHDCoordSys
import MHDRectPM
import MHDSystem
import MHDMagnetometer
import numpy
#import scipy
from scipy.spatial.transform import Rotation as R


class MHDUniverse:

    def __init__(self,dt,tMax):
        self.dt = dt
        self.tMax = tMax
        self.cs = MHDCoordSys.MHDCoordSys()
        self.systems = []
        self.magnets = []

    def addSystem(self,sys):
        self.systems.append(sys)

    def addMagnet(self,mag):
        self.magnets.append(mag)


    def getBFieldForMgtr(self,mgtr):
        #print('getBFieldForMgtr')
        bMgtrSum = numpy.array([0.0, 0.0, 0.0], numpy.double)
        for mag in self.magnets:
            #print('eval mag:%s'%mag)
            sys = mgtr.system
            aMgtr = numpy.array([mgtr.cs.x,mgtr.cs.y,mgtr.cs.z],numpy.double)
            rMgtr = R.from_euler(mgtr.cs.rot, [mgtr.cs.thetaX, mgtr.cs.thetaY, mgtr.cs.thetaZ], degrees=True)
            aSys = rMgtr.apply(aMgtr)
            rSys = R.from_euler(sys.cs.rot, [sys.cs.thetaX, sys.cs.thetaY, sys.cs.thetaZ], degrees=True)
            aUniv = rSys.apply(aSys)
            eSys = numpy.array([sys.cs.x,sys.cs.y,sys.cs.z],numpy.double)
            eUniv = rSys.apply(eSys)
            cUniv = eUniv + aUniv

            dMag = numpy.array([mag.cs.x, mag.cs.y, mag.cs.z], numpy.double)
            rMag = R.from_euler(mag.cs.rot, [mag.cs.thetaX, mag.cs.thetaY, mag.cs.thetaZ], degrees=True)
            dUniv = dMag#rMag.apply(dMag)
            fUniv = dUniv - cUniv

            rInvMag = rMag.inv()
            fMag = rInvMag.apply(fUniv)

            #need to get b field now from fMag
            bMag = mag.getB(fMag)
            bUniv = rMag.apply(bMag)
            rSysInv = rSys.inv()
            bSys = rSysInv.apply(bUniv)
            rMgtrInv = rMgtr.inv()
            bMgtr = rMgtrInv.apply(bSys)
            bMgtrSum+= bMgtr

        return bMgtrSum

    def run(self):
        t = 0
        while t<=self.tMax:
            #update positions
            #update system
            for system in self.systems:
                system.step(t)
            t+=self.dt

if 1 == 1:
    for i1 in range(0,1,1):
        #setup universe
        u = MHDUniverse(1,0)

        #setup magnets
        a = .047
        b = .0215
        h = .0097
        magCS = MHDCoordSys.MHDCoordSys(0.0, 0.0, -0.10, 0, 0, 0.0, u.cs)
        mag = MHDRectPM.MHDRectPM(magCS,a, b, h, 1.25e11*1)
        u.addMagnet(mag)

        #setup system
        csSys1 = MHDCoordSys.MHDCoordSys(.0, .0, .00, .0, 00.0, 0.0, u.cs)
        sys1 = MHDSystem.MHDSystemA(csSys1, u)
        u.addSystem(sys1)

        #setup magnetometers
        mgtr1CS = MHDCoordSys.MHDCoordSys(.0, .0, 0.10, 0.0, 0.0, 0.0, sys1.cs)
        mgtr1 = MHDMagnetometer.MHDMagnetometer(mgtr1CS, sys1)
        #sys1.magnetometers.append(mgtr1)

        mgtr2CS = MHDCoordSys.MHDCoordSys(0.0200, .0, -.0, 0.0, 0.0, 0.0, sys1.cs)
        mgtr2 = MHDMagnetometer.MHDMagnetometer(mgtr2CS, sys1)
        #sys1.magnetometers.append(mgtr2)




        #test magnet
        if 1==1:
            for i in range(-10,10,1):
                testMagCS = MHDCoordSys.MHDCoordSys(0, 0, 0.0, 00.0, 0.0, 0.0, u.cs)
                testMag = MHDRectPM.MHDRectPM(magCS, a, b, h, 1.25e11)
                dx = .025
                pos = [0*i*dx,0*i*dx,0+i*dx]
                print('testMagB:%s pos:%s'%(testMag.getB(pos),pos))

        dx = .2
        var = -1 + dx * i1
        print('var:%s'%var)
        vMagCS = MHDCoordSys.MHDCoordSys(.0, .0, -0.270, 0.0, 0.0, 0.0, sys1.cs)
        vMag = MHDRectPM.MHDRectPM(vMagCS, a, b, h, 1.25e11*1)
        vMag.verbose = 1
        sys1.addMagnet(vMag)

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