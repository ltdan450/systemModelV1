import MHDCoordSys
import MHDMagnetometer
import MHDRectPM
#import MHDUniverse
import math
import numpy
from scipy.spatial.transform import Rotation as R
from scipy import signal
import scipy.optimize as opt
import matplotlib.pyplot as plt
import globalConstants as c
import random


def log(tolog):
    verbose = True
    if verbose:
        print(tolog)

class MHDSystemA:

    def __init__(self,cs,universe):

        self.cs = cs
        self.universe = universe
        self.magnetometers=[]
        self.magnet = 0

        self.mgtr1 = 0
        self.mgtr2 = 0
        self.tMeasLast = -10000.0
        self.measInterval = .010

        self.mgtrPairs = []

        self.xsOut = []
        self.ysOut = []

        self.altCS = MHDCoordSys.MHDCoordSys(0,0,0,0,0,0,self.cs)
        self.posFilt = None
        #log ('need to setup self.magnetModel in MHDSystem')

        self.rXs = []
        self.rYs = []
        self.rZs = []

        self.hpFiltParams = []

    def addMagnet(self,magnet):
        self.magnet = magnet

    def addMagnetometer(self,mgtr):
        self.magnetometers.append(mgtr)
        magEstimatedCS = MHDCoordSys.MHDCoordSys(.1,.1,.1,0,0,0,0)
        mgtr.magnet = MHDRectPM.MHDRectPM(magEstimatedCS, self.magnet.a, self.magnet.b, self.magnet.h, self.magnet.J)
        #mgtr.altCS = self.altCS

        for mgtr1 in self.magnetometers:
            for mgtr2 in self.magnetometers:
                if mgtr1!=mgtr2:
                    addFlag = True
                    for pair in self.mgtrPairs:
                        if mgtr1 in pair and mgtr2 in pair:
                            addFlag = False
                    if addFlag==True:
                        self.mgtrPairs.append([mgtr1,mgtr2])
                        log('added pair %s and %s'%(mgtr1,mgtr2))

    def step(self,t):
        if (t-self.tMeasLast)>self.measInterval:
            self.tMeasLast = t
            self.xsOut.append(t)
            mgtrsSorted = self.getOrderedMgtrList()
            mgtr,b = mgtrsSorted[mgtrsSorted.keys()[0]]
            Rmgtr_mgtrCS = mgtr.magnet.getR(b)
            Rsys_sysCS = [Rmgtr_mgtrCS[0]-mgtr.cs.x,Rmgtr_mgtrCS[1]-mgtr.cs.y,Rmgtr_mgtrCS[2]-mgtr.cs.z]
            rHist = self.updateRHist(Rsys_sysCS)
            results = self.getBreathRate(rHist)
            #print('t:%s res:%s'%(t,results) )
            #self.ysOut.append(b[2])
            #print('rMGTR:%s'%(Rmgtr_mgtrCS))
            return results
        else:
            return 0



    def filtPos(self,R):
        if (self.posFilt != None):
            if self.filtType == 'IIR':
                filtZ, self.zZhp = signal.lfilter(self.posFilt[0],self.posFilt[1],[R[2]],zi=self.zZhp)
                filtZ = filtZ[0]
            else:
                #print('filtZ:')
                filtZ = signal.lfilter(self.posFilt[0], self.posFilt[1], [R[2]])

            if random.random()>.99:
                plt.clf()
                plt.figure(1)
                plt.plot(R[2])
                print(R[2])
                print(R)
                plt.show()


            #if abs(filtZ)>5:
            #    filtZ=0
            return [R[0],R[1],filtZ]
        else:
            return R

    def updateRHist(self,R):
        # want to store pos values but how many?
        # assume 1024 values * 2 bytes per value * 3 dims = 6k

        self.nHist = 1024

        self.rXs.insert(0,R[0])
        self.rYs.insert(0,R[1])
        self.rZs.insert(0,R[2])
        self.rXs = self.rXs[:self.nHist]
        self.rYs = self.rYs[:self.nHist]
        self.rZs = self.rZs[:self.nHist]
        return [self.rXs,self.rYs,self.rZs]

    def getBreathRate(self,bHist):
        #only deal with Z for now, correct later
        xsRev = bHist[2][::-1]
        lHist = len(xsRev)
        method = 'zc-FIRHP'
        if method=='zero crossing' and lHist>50:
            from PyAstronomy import pyaC
            tsRev = []
            for i in range(0,lHist,1):
                tsRev.append(self.measInterval*i)

            #print('tsRev:%s'%tsRev)
            #print('xsRev:%s'%xsRev)

            xc, xi = pyaC.zerocross1d(numpy.array(tsRev),numpy.array(xsRev),getIndices=True)
            #print('xc:%s'%(xc))
            try:
                rate = 60/((xi[-1] - xi[-3])*self.measInterval)
            except:
                rate = 0
            return rate
        if method=='zc-FIRHP':
            yHP = signal.lfilter(self.posFilt[0],self.posFilt[1],xsRev)
            self.ysOut.append(yHP[-1])
            from PyAstronomy import pyaC
            tsRev = []
            for i in range(0,lHist,1):
                tsRev.append(self.measInterval*i)
            xc, xi = pyaC.zerocross1d(numpy.array(tsRev),numpy.array(yHP),getIndices=True)

            try:
                rate = 0.5/((xc[-1] - xc[-2])/2+(xc[-2] - xc[-3])/2)*1.0
                mpi = int(round(1.1 / rate / self.measInterval))

                listTrunc = yHP[-mpi:-1]

                maxVal = max(listTrunc)
                minVal = min(listTrunc)
                amp = (maxVal-minVal)*1.05

            except:
                rate = 0
                amp = 0




            if random.random()>1.997:
                print('bHist:%s'%(bHist))
                print('xsRev:%s' % (xsRev))
                print('yHP:%s' % (yHP))
                #print('xc-1:%s'%xc[-1])

                plt.clf()
                plt.figure(1)
                plt.plot(yHP)
                plt.show()
            #self.ysOut.append(yBP[:-1])
            #self.xsOut.append(len(yBP)*self.measInterval)
            #self.ysOut = yBP

            return [rate,amp]



        elif method=='max-min' and lHist>50:
            print('need to implement max-min method')



    def getRate(self,t,R):
        log('get rate for ')
    def getPosError(self,thetaIn):
        E = [0,0,0]
        if len(thetaIn) == 2:
            thetaInbuff = [thetaIn[0],thetaIn[1],0]
            thetaIn = thetaInbuff
        for pair in self.mgtrPairs:
            mgtr1 = pair[0]
            mgtr2 = pair[1]

            #set up corrective rotations, rotCorr rotates TO magnet CS
            thetaCS = [thetaIn[0],thetaIn[1],thetaIn[2]]
            rot_sys_to_mag = R.from_euler('xyz',thetaCS,degrees=True)
            rot_mag_to_sys  = rot_sys_to_mag.inv()

            #get and correct mag fields
            mgtr1B = mgtr1.measure()
            mgtr1B_magnetCS = rot_sys_to_mag.apply(mgtr1B)
            mgtr2B = mgtr2.measure()
            mgtr2B_magnetCS = rot_sys_to_mag.apply(mgtr2B)

            #get computed positions
            R1_magnetCS = mgtr1.magnet.getR(mgtr1B_magnetCS)
            R2_magnetCS = mgtr2.magnet.getR(mgtr2B_magnetCS)
            R1_sysCS = rot_mag_to_sys.apply(R1_magnetCS)
            R2_sysCS = rot_mag_to_sys.apply(R2_magnetCS)

            log('r1_sys:%s \nr2_sys:%s'%(R1_sysCS,R2_sysCS))

            #get error
            E[0] = E[0]+(mgtr2.cs.x-mgtr1.cs.x)-(R2_sysCS[0]-R1_sysCS[0])
            E[1] = E[1]+(mgtr2.cs.y-mgtr1.cs.y)-(R2_sysCS[1]-R1_sysCS[1])
            E[2] = E[2]+(mgtr2.cs.z-mgtr1.cs.z)-(R2_sysCS[2]-R1_sysCS[2])
            #E = [,,]

            #log('E:%s'%(E))
        eNorm = math.sqrt(E[0]**2+E[1]**2+E[2]**2)
        E = numpy.array(E)
        print('E:%s'%E)
        print('ENorm:%s'%(eNorm))

        return eNorm
        #return E

    def getOrientation(self,guess,tol=1e-6):
        #if guess==0:
        #    guess = [.0,.0,.0]

        #x = opt.anderson(self.getPosError, guess, f_tol=5E-7)
        #x = opt.newton_krylov(self.getPosError, guess, f_tol=tol, x_tol = 1.0, callback=self.itStep)
        x = opt.minimize(self.getPosError,guess[0:3],method = 'Nelder-Mead')

        print('x:%s'%x)

        #x = opt.newton_krylov(self.getPosError, )
        #x = opt.newton(self.getPosError, guess)


        rerun = False
        if 1 ==1:
            for i in range(0,3,1):
                if x[i] > 180:
                    x[i] = 0#math.remainder(x[i],180)*180.0
                    rerun = True

                if x[i] < -180:
                    x[i] = 0#math.remainder(x[i],180)*180.0
                    rerun = True
            if rerun:
                x = self.getOrientation(x)

        print('orientation:%s'%(x))
        return x

    def getOrientationMeth2(self,guess):
        print('need to impelment')

    def getOrderedMgtrList(self):
        method = 'bMagnitude'
        if method == 'bMagnitude':
            magsBuff = {}
            for mgtr in self.magnetometers:
                b = mgtr.measure()
                mag = b[0]**2+b[1]**2+b[2]**2
                magsBuff[mag]=[mgtr, b]
            magsSorted = sorted(magsBuff.keys())
            magsOut = {}
            for mag in magsSorted:
                magsOut[mag]=magsBuff[mag]
            return magsOut


    def itStep(self,x,f):
        print('f:%s x:%s'%(x,f))


    def testEnd(self):
        plt.clf()
        plt.figure(1)
        if 1 == 0:
            bpmThresh = 10.
            fThresh = bpmThresh / 60.0
            fHigh = 2. / 60.
            fLow = 1. / 60.
            ord = 6
            fNyq = 0.5 * 1.0 / (self.measInterval)
            print('fH:%s fL:%s nyq:%s' % (fThresh, fLow, fNyq))
            # b, a = signal.butter(ord,[fLow/fNyq, fHigh/fNyq],btype='band')
            if 1 == 1:
                b, a = signal.butter(ord, [fThresh / fNyq], btype='highpass')
            else:
                b, a = signal.butter(ord, [fHigh / fNyq], btype='highpass')
            yBP = signal.lfilter(b, a, self.ysOut)
            self.ysOut = yBP

        if 1 == 0:
            self.ysOut= self.rZs
        if len(self.ysOut)>0:
            if len(self.ysOut)==len(self.xsOut):
                plt.plot(self.xsOut,self.ysOut)
            else:
                #print('ysOut plot:')
                plt.plot(self.ysOut)
            plt.show()



    def test2(self):
        a = 4*.0254#0.048
        b = 6*.0254#0.022
        h = 0.5*.0254#0.011

        globalCS = MHDCoordSys.MHDCoordSys()

        #setup magnet
        xRm = 0.0
        yRm = 0.0
        zRm = -0.0
        tXRm = 000.00
        tYRm = 010.0
        tZRm = 00.0
        rmCS = MHDCoordSys.MHDCoordSys(xRm, yRm, zRm,tXRm,tYRm,tZRm)
        rm = MHDRectPM.MHDRectPM(rmCS,a,b,h,1.0)
        rm.setJandKForBr(0.31)

        inputs = []
        outputs = []
        out1 = []
        out2 = []
        out3 = []
        #corrective CS
        vMin = -20.0
        vMax = 20.0
        nSamps = 60
        for i in range (0,nSamps,1):
            f = float(i)
            input = float(vMin + f * float((vMax-vMin)/float(nSamps)))
            log('input:%s'%input)
            tXCorr = 0
            tYCorr = input#input/2.0
            tZCorr = 0#input/4.0
            corrCS = MHDCoordSys.MHDCoordSys(0,0,0,tXCorr,tYCorr,tZCorr)
            rotCorr = R.from_euler('xyz',[corrCS.thetaX,corrCS.thetaY,corrCS.thetaZ],degrees=True)
            invRotCorr = rotCorr.inv()


            #setup mgtr1
            x1 =-0.0323
            y1 = 0.0
            z1 = 0.1832
            tX1 = 0.0
            tY1 = 0.0
            tZ1 = 0.0
            cs1 = MHDCoordSys.MHDCoordSys(x1, y1, z1,tX1,tY1,tZ1)
            r1 = [cs1.x-rmCS.x, cs1.y-rmCS.y, cs1.z - rmCS.z]
            r1 = numpy.array(r1)
            log('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            log('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            log('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            log('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorr.apply(b1Real_mgtr1)
            log('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            log('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            log('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))

            #setup mgtr2
            x2 = 0.1315
            y2 = 0.0
            z2 = 0.1315
            tX2 = 0.0
            tY2 = 0.0
            tZ2 = 0.0
            cs2 = MHDCoordSys.MHDCoordSys(x2, y2, z2,tX2,tY2,tZ2)
            r2 = [cs2.x-rmCS.x, cs2.y-rmCS.y, cs2.z - rmCS.z]
            r2 = numpy.array(r2)
            log('r:%s'%r2)
            realDist = (r2[0]**2 + r2[1]**2 + r2[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r2_mag = invRot.apply(r2)
            log('r2_mag:%s'%r2_mag)

            b2Real_mag = rm.getB(r2_mag)
            log('b2Real_mag:%s'%b2Real_mag)
            b2Real_mgtr2 = rot.apply(b2Real_mag)
            log('b2Real_mgtr2:%s'%b2Real_mgtr2)
            b2Real_mgtr2Corr = invRotCorr.apply(b2Real_mgtr2)
            log('b2Real_mgtr2Corr:%s' % b2Real_mgtr2Corr)
            r2Calc_mgtr2 = rm.getR(b2Real_mgtr2Corr)
            log('r2Calc_mgtr1:%s'%r2Calc_mgtr2)
            log('r2Calc_magnitude:%s'%(r2Calc_mgtr2[0]**2 + r2Calc_mgtr2[1]**2 + r2Calc_mgtr2[2]**2)**(0.5))

            #eval positions
            R2R1Real_sys = [cs2.x-cs1.x,cs2.y-cs1.y,cs2.z-cs1.z]
            R2R1RealMag=(R2R1Real_sys[0] ** 2 + R2R1Real_sys[1] ** 2 + R2R1Real_sys[2] ** 2) ** (0.5)
            log('R2R1RealMag:%s' %R2R1RealMag)
            #R2R1Real_sys = rotCorr.apply(R2R1Real_mag)
            log('R2R1Real_sys:%s' % R2R1Real_sys)

            R2R1calc_mag = [r2Calc_mgtr2[0]-rCalc_mgtr1[0],r2Calc_mgtr2[1]-rCalc_mgtr1[1],r2Calc_mgtr2[2]-rCalc_mgtr1[2]]
            R2R1CalcMag = (R2R1calc_mag[0] ** 2 + R2R1calc_mag[1] ** 2 + R2R1calc_mag[2] ** 2) ** (0.5)
            log('R2R1CalcMag:%s' % R2R1CalcMag)
            #might try converting back to global/system coords
            R2R1calc_sys = rotCorr.apply(R2R1calc_mag)
            log('R2R1Calc_sys:%s' % R2R1calc_sys)


            eVec = [R2R1Real_sys[0]-R2R1calc_sys[0],R2R1Real_sys[1]-R2R1calc_sys[1],R2R1Real_sys[2]-R2R1calc_sys[2]]

            error = R2R1RealMag-R2R1CalcMag
            log('error:%s'%error)
            goalError=1.438244776413855e-05*2.0
            errorFactor = error/goalError
            log('errorFactor:%s'%errorFactor)
            inputs.append(input)
            outputs.append(error)
            out1.append(eVec[0])
            out2.append(eVec[1])
            out3.append(eVec[2])

        for i in range(0,len(inputs)):
            log('%s : %s'%(inputs[i],outputs[i]))

        plt.clf()
        plt.figure(1)
        plt.plot(inputs,outputs,label='errorMags')
        plt.plot(inputs, out1, label='errorX')
        plt.plot(inputs, out2, label='errorY')
        plt.plot(inputs, out3, label='errorZ')
        plt.xlabel('Degrees of rotation, y-axis')
        plt.ylabel('Positional error of magnetometers')
        plt.legend()
        plt.grid()
        plt.show()

    def test3(self):
        a = 1.0625*.0254#0.048
        b = 0.75*.0254#0.022
        h = .125*.0254#0.011

        #test set J and K
        #self.setJandKForBr(0.31)

        globalCS = MHDCoordSys.MHDCoordSys()

        #setup magnet
        xRm = 0.0
        yRm = 0.0
        zRm = -0.0
        tXRm = 000.00
        tYRm = 000.0
        tZRm = 00.0
        rmCS = MHDCoordSys.MHDCoordSys(xRm, yRm, zRm,tXRm,tYRm,tZRm)
        rm = MHDRectPM.MHDRectPM(rmCS,a,b,h,1.0)
        rm.setJandKForBr(1.3)

        inputs = []
        outputs = []
        out1 = []
        out2 = []
        out3 = []
        #corrective CS
        vMin =-0.001
        vMax = 0.001
        nSamps = 60
        for i in range (0,nSamps,1):
            f = float(i)
            input = float(vMin + f * float((vMax-vMin)/float(nSamps)))
            log('input:%s'%input)
            tXCorr = 0
            tYCorr = 0#input/2.0
            tZCorr = 0#input/4.0
            corrCS = MHDCoordSys.MHDCoordSys(0,0,0,tXCorr,tYCorr,tZCorr)
            rotCorr = R.from_euler('xyz',[corrCS.thetaX,corrCS.thetaY,corrCS.thetaZ],degrees=True)
            invRotCorr = rotCorr.inv()


            #setup mgtr1
            x1 = 0.0
            y1 = -0.0
            z1 = input + .3
            tX1 = 0.0
            tY1 = 0.0
            tZ1 = 0.0
            cs1 = MHDCoordSys.MHDCoordSys(x1, y1, z1,tX1,tY1,tZ1)
            r1 = [cs1.x-rmCS.x, cs1.y-rmCS.y, cs1.z - rmCS.z]
            r1 = numpy.array(r1)
            log('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            log('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            log('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            log('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorr.apply(b1Real_mgtr1)
            log('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            log('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            log('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))


            inputs.append(input)
            outputs.append(math.sqrt(b1Real_mgtr1[0]**2+b1Real_mgtr1[1]**2+b1Real_mgtr1[2]**2))
            out1.append(b1Real_mgtr1[0]*1e6)
            out2.append(b1Real_mgtr1[1]*1e6)
            out3.append(b1Real_mgtr1[2]*1e6)

        for i in range(0,len(inputs)):
            log('%s : %s'%(inputs[i],outputs[i]))

        plt.clf()
        plt.figure(1)
        #plt.plot(inputs,outputs,label='magMagnitude')
        #plt.plot(inputs, out1, label='magX')
        #plt.plot(inputs, out2, label='magY')
        plt.plot(inputs, out3, label='magZ')
        plt.xlabel('xPos (mm)')
        plt.ylabel('BMag')
        plt.legend()
        plt.grid()
        plt.show()

    def test4(self):
        import MHDUniverse
        #this is an extension of test 2
        a = c.a#4*.0254#0.048
        b = c.b#6*.0254#0.022
        h = c.h#0.5*.0254#0.011

        globalCS = MHDCoordSys.MHDCoordSys()

        univ = MHDUniverse.MHDUniverse(1,1)
        #univ = MHD MHDUniverse(1,1)



        #setup magnet
        xRm = c.xRm#0.50
        yRm = c.yRm#0.20
        zRm = c.zRm#-0.50
        tXRm = c.tXRm#0000.00
        tYRm = c.tYRm#0010.0
        tZRm = c.tZRm#0000.0
        rmCS = MHDCoordSys.MHDCoordSys(xRm, yRm, zRm,tXRm,tYRm,tZRm,univ.cs)
        rm = MHDRectPM.MHDRectPM(rmCS,a,b,h,1.0)
        rm.setJandKForBr(c.Br)

        inputs = []
        outputs = []
        out1 = []
        out2 = []
        out3 = []
        #corrective CS
        vMin = c.tYMin
        vMax = c.tYMax
        nSamps = 2
        for i in range (0,nSamps,1):

            #step through a number of angles and use them to make corrective coordinate systems
            f = float(i)
            input = float(vMin + f * float((vMax-vMin)/float(nSamps)))
            log('input:%s'%input)
            tXCorrective = 0
            tYCorrective = input
            tZCorrective = 0
            correctiveCS = MHDCoordSys.MHDCoordSys(0,0,0,tXCorrective,tYCorrective,tZCorrective)

            #generate corrective rotation and inverse
            rotCorrective = R.from_euler('xyz',[correctiveCS.thetaX,correctiveCS.thetaY,correctiveCS.thetaZ],degrees=True)
            invRotCorrective = rotCorrective.inv()

            csSys1 = MHDCoordSys.MHDCoordSys(0,0,0.1,0,0,0,univ.cs)



            #setup mgtr1
            x1 =c.x1#-0.0323
            y1 = c.y1#0.0
            z1 = c.z1#0.1832
            tX1 = c.tX1#0.0
            tY1 = c.tY1#0.0
            tZ1 = c.tZ1#0.0
            cs1 = MHDCoordSys.MHDCoordSys(x1, y1, z1,tX1,tY1,tZ1,csSys1)
            r1 = [cs1.x-rmCS.x, cs1.y-rmCS.y, cs1.z - rmCS.z]
            r1 = numpy.array(r1)
            log('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            log('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            log('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            log('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorrective.apply(b1Real_mgtr1)
            log('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            log('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            log('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))

            #setup mgtr2
            x2 = c.x2#0.1315
            y2 = c.y2#0.0
            z2 = c.z2#0.1315
            tX2 = c.tX2#0.0
            tY2 = c.tY2#0.0
            tZ2 = c.tZ2#0.0
            cs2 = MHDCoordSys.MHDCoordSys(x2, y2, z2,tX2,tY2,tZ2)
            r2 = [cs2.x-rmCS.x, cs2.y-rmCS.y, cs2.z - rmCS.z]
            r2 = numpy.array(r2)
            log('r:%s'%r2)
            realDist = (r2[0]**2 + r2[1]**2 + r2[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r2_mag = invRot.apply(r2)
            log('r2_mag:%s'%r2_mag)

            b2Real_mag = rm.getB(r2_mag)
            log('b2Real_mag:%s'%b2Real_mag)
            b2Real_mgtr2 = rot.apply(b2Real_mag)
            log('b2Real_mgtr2:%s'%b2Real_mgtr2)
            b2Real_mgtr2Corr = invRotCorrective.apply(b2Real_mgtr2)
            log('b2Real_mgtr2Corr:%s' % b2Real_mgtr2Corr)
            r2Calc_mgtr2 = rm.getR(b2Real_mgtr2Corr)
            log('r2Calc_mgtr1:%s'%r2Calc_mgtr2)
            log('r2Calc_magnitude:%s'%(r2Calc_mgtr2[0]**2 + r2Calc_mgtr2[1]**2 + r2Calc_mgtr2[2]**2)**(0.5))

            #eval positions
            R2R1Real_sys = [cs2.x-cs1.x,cs2.y-cs1.y,cs2.z-cs1.z]
            R2R1RealMag=(R2R1Real_sys[0] ** 2 + R2R1Real_sys[1] ** 2 + R2R1Real_sys[2] ** 2) ** (0.5)
            log('R2R1RealMag:%s' %R2R1RealMag)
            #R2R1Real_sys = rotCorr.apply(R2R1Real_mag)
            log('R2R1Real_sys:%s' % R2R1Real_sys)

            R2R1calc_mag = [r2Calc_mgtr2[0]-rCalc_mgtr1[0],r2Calc_mgtr2[1]-rCalc_mgtr1[1],r2Calc_mgtr2[2]-rCalc_mgtr1[2]]
            R2R1CalcMag = (R2R1calc_mag[0] ** 2 + R2R1calc_mag[1] ** 2 + R2R1calc_mag[2] ** 2) ** (0.5)
            log('R2R1CalcMag:%s' % R2R1CalcMag)
            #might try converting back to global/system coords
            R2R1calc_sys = rotCorrective.apply(R2R1calc_mag)
            log('R2R1Calc_sys:%s' % R2R1calc_sys)


            eVec = [R2R1Real_sys[0]-R2R1calc_sys[0],R2R1Real_sys[1]-R2R1calc_sys[1],R2R1Real_sys[2]-R2R1calc_sys[2]]

            error = R2R1RealMag-R2R1CalcMag
            log('error:%s'%error)
            log('error vector:%s'%eVec)
            goalError=1.438244776413855e-05*2.0
            errorFactor = error/goalError
            log('errorFactor:%s'%errorFactor)
            inputs.append(input)
            #outputs.append(error)
            outputs.append((eVec[0]**2+eVec[1]**2+eVec[2]**2)**0.5)
            #outputs.append(eVec[0] + eVec[1] + eVec[2])
            out1.append(eVec[0])
            out2.append(eVec[1])
            out3.append(eVec[2])

        for i in range(0,len(inputs)):
            log('%s : %s'%(inputs[i],outputs[i]))

        plt.clf()
        plt.figure(1)
        plt.plot(inputs,outputs,label='errorMags')
        plt.plot(inputs, out1, label='errorX')
        plt.plot(inputs, out2, label='errorY')
        plt.plot(inputs, out3, label='errorZ')
        plt.xlabel('Degrees of rotation, y-axis')
        plt.ylabel('Positional error of magnetometers')
        plt.legend()
        plt.grid()
        plt.show()

class MHDSystemB:

    def __init__(self, cs, universe):

        self.cs = cs
        self.universe = universe
        self.magnetometers=[]
        self.magnets = []
        self.magAssbys = []

        self.posComp = numpy.array([1.0,1.0,1.0],numpy.double)
        self.orientationComp = numpy.array([0.0,0.0,0.0],numpy.double)

        self.mgtr1 = 0
        self.mgtr2 = 0
        self.tMeasLast = -10000.0
        self.measInterval = .010

        self.mgtrPairs = []

        self.xsOut = []
        self.ysOut = []

        self.altCS = MHDCoordSys.MHDCoordSys(0,0,0,0,0,0,self.cs)
        self.posFilt = None
        #log ('need to setup self.magnetModel in MHDSystem')

        self.rXs = []
        self.rYs = []
        self.rZs = []

        self.hpFiltParams = []

    def addMagnet(self, magnet):
        self.magnets.append(magnet)
    def addMagAssby(self,magAssby):
        self.magAssbys.append(magAssby)



    def addMagnetometer(self,mgtr):
        self.magnetometers.append(mgtr)

        mgtr.positionComp = numpy.array([1.0, 1.0, 1.0], numpy.double)
        mgtr.orientationComp = numpy.array([0.0, 0.0, 0.0], numpy.double)

        for mgtr1 in self.magnetometers:
            for mgtr2 in self.magnetometers:
                if mgtr1!=mgtr2:
                    addFlag = True
                    for pair in self.mgtrPairs:
                        if mgtr1 in pair and mgtr2 in pair:
                            addFlag = False
                    if addFlag==True:
                        self.mgtrPairs.append([mgtr1,mgtr2])
                        log('added pair %s and %s'%(mgtr1,mgtr2))

    def getMgtrBMeasurement(self,mgtr):
        a_Sys = numpy.array([mgtr.cs.x, mgtr.cs.y, mgtr.cs.z], numpy.double)
        b_Sys = self.universe.getBFieldForA(a_Sys)
        # do something with b_Sys
        rotMgtr = R.from_euler(mgtr.cs.rot, [mgtr.cs.thetaX, mgtr.cs.thetaY, mgtr.cs.thetaZ], degrees=True)
        rInvMgtr = rotMgtr.inv()
        b_mgtr = rInvMgtr.apply(b_Sys)
        return b_mgtr



    def step(self,t):
        if (t-self.tMeasLast)>self.measInterval:
            self.tMeasLast = t
            self.xsOut.append(t)
            mgtrsSorted = self.getOrderedMgtrList()
            #mgtr, b = mgtrsSorted[mgtrsSorted.keys()[0]]
            for mgtr in self.magnetometers:
                b_mgtr = self.getMgtrBMeasurement(mgtr)
                # now solve for R!
                vMagAssby = self.magAssbys[0]
                a_Sys = numpy.array([mgtr.cs.x, mgtr.cs.y, mgtr.cs.z], numpy.double)

                #k_sys needs to be solved/updated
                k_sys = self.posComp
                i_sys = k_sys - a_Sys
                rotMagAssby = R.from_euler(sys.cs.rot, [self.orientationComp[0], self.orientationComp[1], self.orientationComp[2]], degrees=True)
                rInvMagAssby = rotMagAssby.inv()
                i_magAssby = rInvMagAssby.apply(i_sys)
                b_MagAssby = vMagAssby.getBForI(i_magAssby)
                b_sys = rotMagAssby.apply(b_MagAssby)
                rotMgtr = R.from_euler(mgtr.cs.rot, [mgtr.cs.thetaX, mgtr.cs.thetaY, mgtr.cs.thetaZ], degrees=True)
                rInvMgtr = rotMgtr.inv()
                b_mgtr_calc = rInvMgtr.apply(b_sys)
                err = b_mgtr - b_mgtr_calc
                #maybe turn this into a funciton with inputs of k_sys, rotMagAssem, b_meas_mgtr???





            return results
        else:
            return 0

    def updatePosition(self,b_mgtr,guessX,guessY,guessZ):
        guessXP = self.getXP(guessX, guessY, guessZ)
        guessYP = self.getYP(guessX, guessY, guessZ)
        guessZP = self.getZP(guessX, guessY, guessZ)
        posPrime = numpy.array([[guessXP], [guessYP], [guessZP]], numpy.double)

        #Bx *= -1.0
        #By *= -1.0


        BxPG = self.BxP(guessXP, guessYP, guessZP)
        ByPG = self.ByP(guessXP, guessYP, guessZP)
        BzPG = self.BzP(guessXP, guessYP, guessZP)

        log('xPG:%s yPG:%s zPG:%s'%(guessXP,guessYP,guessZP))

        log('BxPG:%s ByPG:%s BzPG:%s'%(BxPG,ByPG,BzPG))
        log('BxP:%s ByP:%s BzP:%s' % (Bx, By, Bz))

        objFunction = numpy.array([[BxPG-Bx],[ByPG-By],[BzPG-Bz]],numpy.double)

        BPErrorNorm = numpy.linalg.norm(objFunction)
        log('BPErrorNorm:%s'%BPErrorNorm)

        trialPos = 0
        validPosFlag = 0
        errorThreshold = 1e-9

        log('begin nonlinear equation solve')

        while (trialPos<=10000):
            log('trialPos loop begin:%i'%trialPos)
            iterations = 0
            while iterations<100:
                log('iterations loop begin:%i'%iterations)

                # If pos > 4m away from sensor, guess new poss
                if numpy.linalg.norm(posPrime)>4.0:
                    guessXP = (random.random() - 0.5) * .20
                    guessYP = (random.random() - 0.5) * .20
                    guessZP = (random.random() - 0.0) * .20
                    posPrime = numpy.array([[guessXP], [guessYP], [guessZP]], numpy.double)

                log('posPrime: %s '%str(posPrime))

                jac = self.getUpdatedJacobian(guessXP,guessYP,guessZP,.000001)
                log('jacobian: %s'%str(jac))

                invJac = numpy.linalg.inv(jac)
                log('inverse jac: %s'%str(invJac))

                delta = invJac.dot(objFunction)
                log('delta: %s'%str(delta))

                posPrime = numpy.subtract(posPrime,delta)
                log('updatedPos:%s'%posPrime)

                guessXP = posPrime[0,0]
                guessYP = posPrime[1,0]
                guessZP = posPrime[2,0]

                if guessZP < 0.0:
                    guessZP*=1.0

                BxPG = self.BxP(guessXP, guessYP, guessZP)
                ByPG = self.ByP(guessXP, guessYP, guessZP)
                BzPG = self.BzP(guessXP, guessYP, guessZP)
                log('BxPG:%s ByPG:%s BzPG:%s' % (BxPG, ByPG, BzPG))
                log('BxP:%s ByP:%s BzP:%s' % (Bx, By, Bz))
                objFunction = numpy.array([[BxPG - Bx], [ByPG - By], [BzPG - Bz]], numpy.double)
                BPErrorNorm = numpy.linalg.norm(objFunction)
                log('error norm:%f'%(BPErrorNorm))

                x = posPrime[0, 0]
                y = posPrime[1, 0]
                z = posPrime[2, 0]

                if (BPErrorNorm < errorThreshold):

                    validPosFlag = 1
                    break                                       #break iterations loop
                iterations+=1
            #proceed
            log('validPosFlag:%i'%validPosFlag)
            if self.getZ(guessXP, guessYP, guessZP)<0:
                validPosFlag = 0
            if validPosFlag == 1:
                guessXP = posPrime[0,0]
                guessYP = posPrime[1,0]
                guessZP = posPrime[2,0]

                self.x = self.getX(guessXP, guessYP, guessZP)
                self.y = self.getY(guessXP, guessYP, guessZP)
                self.z = self.getZ(guessXP, guessYP, guessZP)

                self.posVect = numpy.array([[self.x],[self.y],[self.z]],numpy.double)
                break                                           #break trialpos loop

            else:
                trialPos+=1
                guessXP = (random.random() - 0.5) * .20
                guessYP = (random.random() - 0.5) * .20
                guessZP = (random.random() - 0.0) * .20

                posPrime = numpy.array([[guessXP], [guessYP], [guessZP]], numpy.double)

        if validPosFlag == 0:
            log('ERROR: could not resolve positon')
        else:
            return [self.x,self.y,self.z]



    def filtPos(self,R):
        if (self.posFilt != None):
            if self.filtType == 'IIR':
                filtZ, self.zZhp = signal.lfilter(self.posFilt[0],self.posFilt[1],[R[2]],zi=self.zZhp)
                filtZ = filtZ[0]
            else:
                #print('filtZ:')
                filtZ = signal.lfilter(self.posFilt[0], self.posFilt[1], [R[2]])

            if random.random()>.99:
                plt.clf()
                plt.figure(1)
                plt.plot(R[2])
                print(R[2])
                print(R)
                plt.show()


            #if abs(filtZ)>5:
            #    filtZ=0
            return [R[0],R[1],filtZ]
        else:
            return R

    def updateRHist(self,R):
        # want to store pos values but how many?
        # assume 1024 values * 2 bytes per value * 3 dims = 6k

        self.nHist = 1024

        self.rXs.insert(0,R[0])
        self.rYs.insert(0,R[1])
        self.rZs.insert(0,R[2])
        self.rXs = self.rXs[:self.nHist]
        self.rYs = self.rYs[:self.nHist]
        self.rZs = self.rZs[:self.nHist]
        return [self.rXs,self.rYs,self.rZs]

    def getBreathRate(self,bHist):
        #only deal with Z for now, correct later
        xsRev = bHist[2][::-1]
        lHist = len(xsRev)
        method = 'zc-FIRHP'
        if method=='zero crossing' and lHist>50:
            from PyAstronomy import pyaC
            tsRev = []
            for i in range(0,lHist,1):
                tsRev.append(self.measInterval*i)

            #print('tsRev:%s'%tsRev)
            #print('xsRev:%s'%xsRev)

            xc, xi = pyaC.zerocross1d(numpy.array(tsRev),numpy.array(xsRev),getIndices=True)
            #print('xc:%s'%(xc))
            try:
                rate = 60/((xi[-1] - xi[-3])*self.measInterval)
            except:
                rate = 0
            return rate
        if method=='zc-FIRHP':
            yHP = signal.lfilter(self.posFilt[0],self.posFilt[1],xsRev)
            self.ysOut.append(yHP[-1])
            from PyAstronomy import pyaC
            tsRev = []
            for i in range(0,lHist,1):
                tsRev.append(self.measInterval*i)
            xc, xi = pyaC.zerocross1d(numpy.array(tsRev),numpy.array(yHP),getIndices=True)

            try:
                rate = 0.5/((xc[-1] - xc[-2])/2+(xc[-2] - xc[-3])/2)*1.0
                mpi = int(round(1.1 / rate / self.measInterval))

                listTrunc = yHP[-mpi:-1]

                maxVal = max(listTrunc)
                minVal = min(listTrunc)
                amp = (maxVal-minVal)*1.05

            except:
                rate = 0
                amp = 0




            if random.random()>1.997:
                print('bHist:%s'%(bHist))
                print('xsRev:%s' % (xsRev))
                print('yHP:%s' % (yHP))
                #print('xc-1:%s'%xc[-1])

                plt.clf()
                plt.figure(1)
                plt.plot(yHP)
                plt.show()
            #self.ysOut.append(yBP[:-1])
            #self.xsOut.append(len(yBP)*self.measInterval)
            #self.ysOut = yBP

            return [rate,amp]



        elif method=='max-min' and lHist>50:
            print('need to implement max-min method')



    def getRate(self,t,R):
        log('get rate for ')
    def getPosError(self,thetaIn):
        E = [0,0,0]
        if len(thetaIn) == 2:
            thetaInbuff = [thetaIn[0],thetaIn[1],0]
            thetaIn = thetaInbuff
        for pair in self.mgtrPairs:
            mgtr1 = pair[0]
            mgtr2 = pair[1]

            #set up corrective rotations, rotCorr rotates TO magnet CS
            thetaCS = [thetaIn[0],thetaIn[1],thetaIn[2]]
            rot_sys_to_mag = R.from_euler('xyz',thetaCS,degrees=True)
            rot_mag_to_sys  = rot_sys_to_mag.inv()

            #get and correct mag fields
            mgtr1B = mgtr1.measure()
            mgtr1B_magnetCS = rot_sys_to_mag.apply(mgtr1B)
            mgtr2B = mgtr2.measure()
            mgtr2B_magnetCS = rot_sys_to_mag.apply(mgtr2B)

            #get computed positions
            R1_magnetCS = mgtr1.magnet.getR(mgtr1B_magnetCS)
            R2_magnetCS = mgtr2.magnet.getR(mgtr2B_magnetCS)
            R1_sysCS = rot_mag_to_sys.apply(R1_magnetCS)
            R2_sysCS = rot_mag_to_sys.apply(R2_magnetCS)

            log('r1_sys:%s \nr2_sys:%s'%(R1_sysCS,R2_sysCS))

            #get error
            E[0] = E[0]+(mgtr2.cs.x-mgtr1.cs.x)-(R2_sysCS[0]-R1_sysCS[0])
            E[1] = E[1]+(mgtr2.cs.y-mgtr1.cs.y)-(R2_sysCS[1]-R1_sysCS[1])
            E[2] = E[2]+(mgtr2.cs.z-mgtr1.cs.z)-(R2_sysCS[2]-R1_sysCS[2])
            #E = [,,]

            #log('E:%s'%(E))
        eNorm = math.sqrt(E[0]**2+E[1]**2+E[2]**2)
        E = numpy.array(E)
        print('E:%s'%E)
        print('ENorm:%s'%(eNorm))

        return eNorm
        #return E

    def getOrientation(self,guess,tol=1e-6):
        #if guess==0:
        #    guess = [.0,.0,.0]

        #x = opt.anderson(self.getPosError, guess, f_tol=5E-7)
        #x = opt.newton_krylov(self.getPosError, guess, f_tol=tol, x_tol = 1.0, callback=self.itStep)
        x = opt.minimize(self.getPosError,guess[0:3],method = 'Nelder-Mead')

        print('x:%s'%x)

        #x = opt.newton_krylov(self.getPosError, )
        #x = opt.newton(self.getPosError, guess)


        rerun = False
        if 1 ==1:
            for i in range(0,3,1):
                if x[i] > 180:
                    x[i] = 0#math.remainder(x[i],180)*180.0
                    rerun = True

                if x[i] < -180:
                    x[i] = 0#math.remainder(x[i],180)*180.0
                    rerun = True
            if rerun:
                x = self.getOrientation(x)

        print('orientation:%s'%(x))
        return x

    def getOrientationMeth2(self,guess):
        print('need to impelment')

    def getOrderedMgtrList(self):
        method = 'bMagnitude'
        if method == 'bMagnitude':
            magsBuff = {}
            for mgtr in self.magnetometers:
                b = mgtr.measure()
                mag = b[0]**2+b[1]**2+b[2]**2
                magsBuff[mag]=[mgtr, b]
            magsSorted = sorted(magsBuff.keys())
            magsOut = {}
            for mag in magsSorted:
                magsOut[mag]=magsBuff[mag]
            return magsOut


    def itStep(self,x,f):
        print('f:%s x:%s'%(x,f))


    def testEnd(self):
        plt.clf()
        plt.figure(1)
        if 1 == 0:
            bpmThresh = 10.
            fThresh = bpmThresh / 60.0
            fHigh = 2. / 60.
            fLow = 1. / 60.
            ord = 6
            fNyq = 0.5 * 1.0 / (self.measInterval)
            print('fH:%s fL:%s nyq:%s' % (fThresh, fLow, fNyq))
            # b, a = signal.butter(ord,[fLow/fNyq, fHigh/fNyq],btype='band')
            if 1 == 1:
                b, a = signal.butter(ord, [fThresh / fNyq], btype='highpass')
            else:
                b, a = signal.butter(ord, [fHigh / fNyq], btype='highpass')
            yBP = signal.lfilter(b, a, self.ysOut)
            self.ysOut = yBP

        if 1 == 0:
            self.ysOut= self.rZs
        if len(self.ysOut)>0:
            if len(self.ysOut)==len(self.xsOut):
                plt.plot(self.xsOut,self.ysOut)
            else:
                #print('ysOut plot:')
                plt.plot(self.ysOut)
            plt.show()



    def test2(self):
        a = 4*.0254#0.048
        b = 6*.0254#0.022
        h = 0.5*.0254#0.011

        globalCS = MHDCoordSys.MHDCoordSys()

        #setup magnet
        xRm = 0.0
        yRm = 0.0
        zRm = -0.0
        tXRm = 000.00
        tYRm = 010.0
        tZRm = 00.0
        rmCS = MHDCoordSys.MHDCoordSys(xRm, yRm, zRm,tXRm,tYRm,tZRm)
        rm = MHDRectPM.MHDRectPM(rmCS,a,b,h,1.0)
        rm.setJandKForBr(0.31)

        inputs = []
        outputs = []
        out1 = []
        out2 = []
        out3 = []
        #corrective CS
        vMin = -20.0
        vMax = 20.0
        nSamps = 60
        for i in range (0,nSamps,1):
            f = float(i)
            input = float(vMin + f * float((vMax-vMin)/float(nSamps)))
            log('input:%s'%input)
            tXCorr = 0
            tYCorr = input#input/2.0
            tZCorr = 0#input/4.0
            corrCS = MHDCoordSys.MHDCoordSys(0,0,0,tXCorr,tYCorr,tZCorr)
            rotCorr = R.from_euler('xyz',[corrCS.thetaX,corrCS.thetaY,corrCS.thetaZ],degrees=True)
            invRotCorr = rotCorr.inv()


            #setup mgtr1
            x1 =-0.0323
            y1 = 0.0
            z1 = 0.1832
            tX1 = 0.0
            tY1 = 0.0
            tZ1 = 0.0
            cs1 = MHDCoordSys.MHDCoordSys(x1, y1, z1,tX1,tY1,tZ1)
            r1 = [cs1.x-rmCS.x, cs1.y-rmCS.y, cs1.z - rmCS.z]
            r1 = numpy.array(r1)
            log('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            log('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            log('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            log('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorr.apply(b1Real_mgtr1)
            log('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            log('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            log('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))

            #setup mgtr2
            x2 = 0.1315
            y2 = 0.0
            z2 = 0.1315
            tX2 = 0.0
            tY2 = 0.0
            tZ2 = 0.0
            cs2 = MHDCoordSys.MHDCoordSys(x2, y2, z2,tX2,tY2,tZ2)
            r2 = [cs2.x-rmCS.x, cs2.y-rmCS.y, cs2.z - rmCS.z]
            r2 = numpy.array(r2)
            log('r:%s'%r2)
            realDist = (r2[0]**2 + r2[1]**2 + r2[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r2_mag = invRot.apply(r2)
            log('r2_mag:%s'%r2_mag)

            b2Real_mag = rm.getB(r2_mag)
            log('b2Real_mag:%s'%b2Real_mag)
            b2Real_mgtr2 = rot.apply(b2Real_mag)
            log('b2Real_mgtr2:%s'%b2Real_mgtr2)
            b2Real_mgtr2Corr = invRotCorr.apply(b2Real_mgtr2)
            log('b2Real_mgtr2Corr:%s' % b2Real_mgtr2Corr)
            r2Calc_mgtr2 = rm.getR(b2Real_mgtr2Corr)
            log('r2Calc_mgtr1:%s'%r2Calc_mgtr2)
            log('r2Calc_magnitude:%s'%(r2Calc_mgtr2[0]**2 + r2Calc_mgtr2[1]**2 + r2Calc_mgtr2[2]**2)**(0.5))

            #eval positions
            R2R1Real_sys = [cs2.x-cs1.x,cs2.y-cs1.y,cs2.z-cs1.z]
            R2R1RealMag=(R2R1Real_sys[0] ** 2 + R2R1Real_sys[1] ** 2 + R2R1Real_sys[2] ** 2) ** (0.5)
            log('R2R1RealMag:%s' %R2R1RealMag)
            #R2R1Real_sys = rotCorr.apply(R2R1Real_mag)
            log('R2R1Real_sys:%s' % R2R1Real_sys)

            R2R1calc_mag = [r2Calc_mgtr2[0]-rCalc_mgtr1[0],r2Calc_mgtr2[1]-rCalc_mgtr1[1],r2Calc_mgtr2[2]-rCalc_mgtr1[2]]
            R2R1CalcMag = (R2R1calc_mag[0] ** 2 + R2R1calc_mag[1] ** 2 + R2R1calc_mag[2] ** 2) ** (0.5)
            log('R2R1CalcMag:%s' % R2R1CalcMag)
            #might try converting back to global/system coords
            R2R1calc_sys = rotCorr.apply(R2R1calc_mag)
            log('R2R1Calc_sys:%s' % R2R1calc_sys)


            eVec = [R2R1Real_sys[0]-R2R1calc_sys[0],R2R1Real_sys[1]-R2R1calc_sys[1],R2R1Real_sys[2]-R2R1calc_sys[2]]

            error = R2R1RealMag-R2R1CalcMag
            log('error:%s'%error)
            goalError=1.438244776413855e-05*2.0
            errorFactor = error/goalError
            log('errorFactor:%s'%errorFactor)
            inputs.append(input)
            outputs.append(error)
            out1.append(eVec[0])
            out2.append(eVec[1])
            out3.append(eVec[2])

        for i in range(0,len(inputs)):
            log('%s : %s'%(inputs[i],outputs[i]))

        plt.clf()
        plt.figure(1)
        plt.plot(inputs,outputs,label='errorMags')
        plt.plot(inputs, out1, label='errorX')
        plt.plot(inputs, out2, label='errorY')
        plt.plot(inputs, out3, label='errorZ')
        plt.xlabel('Degrees of rotation, y-axis')
        plt.ylabel('Positional error of magnetometers')
        plt.legend()
        plt.grid()
        plt.show()

    def test3(self):
        a = 1.0625*.0254#0.048
        b = 0.75*.0254#0.022
        h = .125*.0254#0.011

        #test set J and K
        #self.setJandKForBr(0.31)

        globalCS = MHDCoordSys.MHDCoordSys()

        #setup magnet
        xRm = 0.0
        yRm = 0.0
        zRm = -0.0
        tXRm = 000.00
        tYRm = 000.0
        tZRm = 00.0
        rmCS = MHDCoordSys.MHDCoordSys(xRm, yRm, zRm,tXRm,tYRm,tZRm)
        rm = MHDRectPM.MHDRectPM(rmCS,a,b,h,1.0)
        rm.setJandKForBr(1.3)

        inputs = []
        outputs = []
        out1 = []
        out2 = []
        out3 = []
        #corrective CS
        vMin =-0.001
        vMax = 0.001
        nSamps = 60
        for i in range (0,nSamps,1):
            f = float(i)
            input = float(vMin + f * float((vMax-vMin)/float(nSamps)))
            log('input:%s'%input)
            tXCorr = 0
            tYCorr = 0#input/2.0
            tZCorr = 0#input/4.0
            corrCS = MHDCoordSys.MHDCoordSys(0,0,0,tXCorr,tYCorr,tZCorr)
            rotCorr = R.from_euler('xyz',[corrCS.thetaX,corrCS.thetaY,corrCS.thetaZ],degrees=True)
            invRotCorr = rotCorr.inv()


            #setup mgtr1
            x1 = 0.0
            y1 = -0.0
            z1 = input + .3
            tX1 = 0.0
            tY1 = 0.0
            tZ1 = 0.0
            cs1 = MHDCoordSys.MHDCoordSys(x1, y1, z1,tX1,tY1,tZ1)
            r1 = [cs1.x-rmCS.x, cs1.y-rmCS.y, cs1.z - rmCS.z]
            r1 = numpy.array(r1)
            log('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            log('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            log('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            log('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorr.apply(b1Real_mgtr1)
            log('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            log('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            log('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))


            inputs.append(input)
            outputs.append(math.sqrt(b1Real_mgtr1[0]**2+b1Real_mgtr1[1]**2+b1Real_mgtr1[2]**2))
            out1.append(b1Real_mgtr1[0]*1e6)
            out2.append(b1Real_mgtr1[1]*1e6)
            out3.append(b1Real_mgtr1[2]*1e6)

        for i in range(0,len(inputs)):
            log('%s : %s'%(inputs[i],outputs[i]))

        plt.clf()
        plt.figure(1)
        #plt.plot(inputs,outputs,label='magMagnitude')
        #plt.plot(inputs, out1, label='magX')
        #plt.plot(inputs, out2, label='magY')
        plt.plot(inputs, out3, label='magZ')
        plt.xlabel('xPos (mm)')
        plt.ylabel('BMag')
        plt.legend()
        plt.grid()
        plt.show()

    def test4(self):
        import MHDUniverse
        #this is an extension of test 2
        a = c.a#4*.0254#0.048
        b = c.b#6*.0254#0.022
        h = c.h#0.5*.0254#0.011

        globalCS = MHDCoordSys.MHDCoordSys()

        univ = MHDUniverse.MHDUniverse(1,1)
        #univ = MHD MHDUniverse(1,1)



        #setup magnet
        xRm = c.xRm#0.50
        yRm = c.yRm#0.20
        zRm = c.zRm#-0.50
        tXRm = c.tXRm#0000.00
        tYRm = c.tYRm#0010.0
        tZRm = c.tZRm#0000.0
        rmCS = MHDCoordSys.MHDCoordSys(xRm, yRm, zRm,tXRm,tYRm,tZRm,univ.cs)
        rm = MHDRectPM.MHDRectPM(rmCS,a,b,h,1.0)
        rm.setJandKForBr(c.Br)

        inputs = []
        outputs = []
        out1 = []
        out2 = []
        out3 = []
        #corrective CS
        vMin = c.tYMin
        vMax = c.tYMax
        nSamps = 2
        for i in range (0,nSamps,1):

            #step through a number of angles and use them to make corrective coordinate systems
            f = float(i)
            input = float(vMin + f * float((vMax-vMin)/float(nSamps)))
            log('input:%s'%input)
            tXCorrective = 0
            tYCorrective = input
            tZCorrective = 0
            correctiveCS = MHDCoordSys.MHDCoordSys(0,0,0,tXCorrective,tYCorrective,tZCorrective)

            #generate corrective rotation and inverse
            rotCorrective = R.from_euler('xyz',[correctiveCS.thetaX,correctiveCS.thetaY,correctiveCS.thetaZ],degrees=True)
            invRotCorrective = rotCorrective.inv()

            csSys1 = MHDCoordSys.MHDCoordSys(0,0,0.1,0,0,0,univ.cs)



            #setup mgtr1
            x1 =c.x1#-0.0323
            y1 = c.y1#0.0
            z1 = c.z1#0.1832
            tX1 = c.tX1#0.0
            tY1 = c.tY1#0.0
            tZ1 = c.tZ1#0.0
            cs1 = MHDCoordSys.MHDCoordSys(x1, y1, z1,tX1,tY1,tZ1,csSys1)
            r1 = [cs1.x-rmCS.x, cs1.y-rmCS.y, cs1.z - rmCS.z]
            r1 = numpy.array(r1)
            log('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            log('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            log('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            log('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorrective.apply(b1Real_mgtr1)
            log('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            log('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            log('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))

            #setup mgtr2
            x2 = c.x2#0.1315
            y2 = c.y2#0.0
            z2 = c.z2#0.1315
            tX2 = c.tX2#0.0
            tY2 = c.tY2#0.0
            tZ2 = c.tZ2#0.0
            cs2 = MHDCoordSys.MHDCoordSys(x2, y2, z2,tX2,tY2,tZ2)
            r2 = [cs2.x-rmCS.x, cs2.y-rmCS.y, cs2.z - rmCS.z]
            r2 = numpy.array(r2)
            log('r:%s'%r2)
            realDist = (r2[0]**2 + r2[1]**2 + r2[2]**2)**(0.5)
            log('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r2_mag = invRot.apply(r2)
            log('r2_mag:%s'%r2_mag)

            b2Real_mag = rm.getB(r2_mag)
            log('b2Real_mag:%s'%b2Real_mag)
            b2Real_mgtr2 = rot.apply(b2Real_mag)
            log('b2Real_mgtr2:%s'%b2Real_mgtr2)
            b2Real_mgtr2Corr = invRotCorrective.apply(b2Real_mgtr2)
            log('b2Real_mgtr2Corr:%s' % b2Real_mgtr2Corr)
            r2Calc_mgtr2 = rm.getR(b2Real_mgtr2Corr)
            log('r2Calc_mgtr1:%s'%r2Calc_mgtr2)
            log('r2Calc_magnitude:%s'%(r2Calc_mgtr2[0]**2 + r2Calc_mgtr2[1]**2 + r2Calc_mgtr2[2]**2)**(0.5))

            #eval positions
            R2R1Real_sys = [cs2.x-cs1.x,cs2.y-cs1.y,cs2.z-cs1.z]
            R2R1RealMag=(R2R1Real_sys[0] ** 2 + R2R1Real_sys[1] ** 2 + R2R1Real_sys[2] ** 2) ** (0.5)
            log('R2R1RealMag:%s' %R2R1RealMag)
            #R2R1Real_sys = rotCorr.apply(R2R1Real_mag)
            log('R2R1Real_sys:%s' % R2R1Real_sys)

            R2R1calc_mag = [r2Calc_mgtr2[0]-rCalc_mgtr1[0],r2Calc_mgtr2[1]-rCalc_mgtr1[1],r2Calc_mgtr2[2]-rCalc_mgtr1[2]]
            R2R1CalcMag = (R2R1calc_mag[0] ** 2 + R2R1calc_mag[1] ** 2 + R2R1calc_mag[2] ** 2) ** (0.5)
            log('R2R1CalcMag:%s' % R2R1CalcMag)
            #might try converting back to global/system coords
            R2R1calc_sys = rotCorrective.apply(R2R1calc_mag)
            log('R2R1Calc_sys:%s' % R2R1calc_sys)


            eVec = [R2R1Real_sys[0]-R2R1calc_sys[0],R2R1Real_sys[1]-R2R1calc_sys[1],R2R1Real_sys[2]-R2R1calc_sys[2]]

            error = R2R1RealMag-R2R1CalcMag
            log('error:%s'%error)
            log('error vector:%s'%eVec)
            goalError=1.438244776413855e-05*2.0
            errorFactor = error/goalError
            log('errorFactor:%s'%errorFactor)
            inputs.append(input)
            #outputs.append(error)
            outputs.append((eVec[0]**2+eVec[1]**2+eVec[2]**2)**0.5)
            #outputs.append(eVec[0] + eVec[1] + eVec[2])
            out1.append(eVec[0])
            out2.append(eVec[1])
            out3.append(eVec[2])

        for i in range(0,len(inputs)):
            log('%s : %s'%(inputs[i],outputs[i]))

        plt.clf()
        plt.figure(1)
        plt.plot(inputs,outputs,label='errorMags')
        plt.plot(inputs, out1, label='errorX')
        plt.plot(inputs, out2, label='errorY')
        plt.plot(inputs, out3, label='errorZ')
        plt.xlabel('Degrees of rotation, y-axis')
        plt.ylabel('Positional error of magnetometers')
        plt.legend()
        plt.grid()
        plt.show()






tests = c.tests

if 2 in tests:
    sys = MHDSystemA(0,0)
    sys.test2()

if 3 in tests:
    sys = MHDSystemA(0,0)
    sys.test3()

if 4 in tests:
    sys = MHDSystemA(0,0)
    sys.test4()

if 5 in tests:
    fSample = 100.0
    tMeas = 1.0/fSample

    fSig = 1.
    dc = .3
    ac = .002

    time = numpy.arange(0, 30, tMeas);
    sig = numpy.sin(2.0*fSig*math.pi*(time))
    sig *= ac
    sig += dc

    plt.clf()
    plt.figure(1)

    while 1:
        plt.clf()
        plt.figure(1)
        paramsIn = raw_input('Enter bpmThresh(,N=6):').split(',')
        print('paramsIn:%s'%paramsIn)
        if len(paramsIn)>1:
            bpm = float(paramsIn[0])
            n = int(paramsIn[1])
        else:
            bpm = float(paramsIn[0])
            n = 5
        bpmThresh = bpm
        fThresh = bpmThresh / 60.0
        fHigh = 2. / 60.
        fLow = 1. / 60.
        ord = n
        fNyq = 0.5 * 1.0 / tMeas
        print('fH:%s fL:%s nyq:%s' % (fThresh, fLow, fNyq))
        # b, a = signal.butter(ord,[fLow/fNyq, fHigh/fNyq],btype='band')
        if 1 == 1:
            b, a = signal.butter(ord, [fThresh / fNyq], btype='highpass')

        yBP = signal.lfilter(b, a, sig)


        if len(time) == len(yBP):
            plt.plot(time, yBP)
        else:
            plt.plot(yBP)

        plt.show()













