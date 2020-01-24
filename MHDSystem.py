import MHDCoordSys
import MHDMagnetometer
import MHDRectPM
#import MHDUniverse
import math
import numpy
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt

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

        xs = []
        ys = []

        varMin = 0
        varMax = varMin * 1 + 0
        dVar = 1
        var = varMin
        while var<=varMax:

            #print('o iterations:%s'%oIterations)
            modCS = MHDCoordSys.MHDCoordSys(0,0,0,var*1,var*1,var*1,self.cs)
            for mgtr in self.magnetometers:
                mgtr.altCS = modCS
                bMeasMgtr = mgtr.measure()      #what mgtr measures in its own coordinate system
                rMgtr = R.from_euler(mgtr.cs.rot, [mgtr.cs.thetaX, mgtr.cs.thetaY, mgtr.cs.thetaZ], degrees=True)
                bMeasSys = rMgtr.apply(bMeasMgtr)
                rMag = R.from_euler(self.magnet.cs.rot, [self.magnet.cs.thetaX, self.magnet.cs.thetaY, self.magnet.cs.thetaZ], degrees=True)
                rMagInv = rMag.inv()
                bMeasMag = rMagInv.apply(bMeasSys)
                #print('bMeasMag:%s'%bMeasMag)
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

                        #print('actualDistance: %s'%actualDistance)
                        #print('measuredDistance: %s'%measuredDistance)
                        #print('errorDistance: %s'%errorDistance)
                        orientationError+=math.sqrt(errorDistance[0]**2+errorDistance[1]**2+errorDistance[2]**2)

            print ('orientationError:%s'%orientationError)
            output = orientationError
            xs.append(var)
            ys.append(output)
            var += dVar

        #print(xs)
        #print(ys)
        plt.figure(1)
        plt.clf()
        plt.plot(xs,ys)
        #plt.show()





            #print('Bmeas:%s'%Bmeas)

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
            print('input:%s'%input)
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
            print('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            print('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            print('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            print('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            print('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorr.apply(b1Real_mgtr1)
            print('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            print('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            print('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))

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
            print('r:%s'%r2)
            realDist = (r2[0]**2 + r2[1]**2 + r2[2]**2)**(0.5)
            print('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r2_mag = invRot.apply(r2)
            print('r2_mag:%s'%r2_mag)

            b2Real_mag = rm.getB(r2_mag)
            print('b2Real_mag:%s'%b2Real_mag)
            b2Real_mgtr2 = rot.apply(b2Real_mag)
            print('b2Real_mgtr2:%s'%b2Real_mgtr2)
            b2Real_mgtr2Corr = invRotCorr.apply(b2Real_mgtr2)
            print('b2Real_mgtr2Corr:%s' % b2Real_mgtr2Corr)
            r2Calc_mgtr2 = rm.getR(b2Real_mgtr2Corr)
            print('r2Calc_mgtr1:%s'%r2Calc_mgtr2)
            print('r2Calc_magnitude:%s'%(r2Calc_mgtr2[0]**2 + r2Calc_mgtr2[1]**2 + r2Calc_mgtr2[2]**2)**(0.5))

            #eval positions
            R2R1Real_sys = [cs2.x-cs1.x,cs2.y-cs1.y,cs2.z-cs1.z]
            R2R1RealMag=(R2R1Real_sys[0] ** 2 + R2R1Real_sys[1] ** 2 + R2R1Real_sys[2] ** 2) ** (0.5)
            print('R2R1RealMag:%s' %R2R1RealMag)
            #R2R1Real_sys = rotCorr.apply(R2R1Real_mag)
            print('R2R1Real_sys:%s' % R2R1Real_sys)

            R2R1calc_mag = [r2Calc_mgtr2[0]-rCalc_mgtr1[0],r2Calc_mgtr2[1]-rCalc_mgtr1[1],r2Calc_mgtr2[2]-rCalc_mgtr1[2]]
            R2R1CalcMag = (R2R1calc_mag[0] ** 2 + R2R1calc_mag[1] ** 2 + R2R1calc_mag[2] ** 2) ** (0.5)
            print('R2R1CalcMag:%s' % R2R1CalcMag)
            #might try converting back to global/system coords
            R2R1calc_sys = rotCorr.apply(R2R1calc_mag)
            print('R2R1Calc_sys:%s' % R2R1calc_sys)


            eVec = [R2R1Real_sys[0]-R2R1calc_sys[0],R2R1Real_sys[1]-R2R1calc_sys[1],R2R1Real_sys[2]-R2R1calc_sys[2]]

            error = R2R1RealMag-R2R1CalcMag
            print('error:%s'%error)
            goalError=1.438244776413855e-05*2.0
            errorFactor = error/goalError
            print('errorFactor:%s'%errorFactor)
            inputs.append(input)
            outputs.append(error)
            out1.append(eVec[0])
            out2.append(eVec[1])
            out3.append(eVec[2])

        for i in range(0,len(inputs)):
            print('%s : %s'%(inputs[i],outputs[i]))

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
            print('input:%s'%input)
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
            print('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            print('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            print('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            print('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            print('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorr.apply(b1Real_mgtr1)
            print('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            print('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            print('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))


            inputs.append(input)
            outputs.append(math.sqrt(b1Real_mgtr1[0]**2+b1Real_mgtr1[1]**2+b1Real_mgtr1[2]**2))
            out1.append(b1Real_mgtr1[0]*1e6)
            out2.append(b1Real_mgtr1[1]*1e6)
            out3.append(b1Real_mgtr1[2]*1e6)

        for i in range(0,len(inputs)):
            print('%s : %s'%(inputs[i],outputs[i]))

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
        a = 4*.0254#0.048
        b = 6*.0254#0.022
        h = 0.5*.0254#0.011

        globalCS = MHDCoordSys.MHDCoordSys()
        univ = MHDUniverse.MHDUniverse(1,1)


        #setup magnet
        xRm = 0.0
        yRm = 0.0
        zRm = -0.10
        tXRm = 000.00
        tYRm = 010.0
        tZRm = 00.0
        rmCS = MHDCoordSys.MHDCoordSys(xRm, yRm, zRm,tXRm,tYRm,tZRm,univ.cs)
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
            print('input:%s'%input)
            tXCorrective = 0
            tYCorrective = input
            tZCorrective = 0
            correctiveCS = MHDCoordSys.MHDCoordSys(0,0,0,tXCorrective,tYCorrective,tZCorrective)
            rotCorrective = R.from_euler('xyz',[correctiveCS.thetaX,correctiveCS.thetaY,correctiveCS.thetaZ],degrees=True)
            invRotCorrective = rotCorrective.inv()

            csSys1 = MHDCoordSys.MHDCoordSys(0,0,0.1,0,0,0,univ.cs)



            #setup mgtr1
            x1 =-0.0323
            y1 = 0.0
            z1 = 0.1832
            tX1 = 0.0
            tY1 = 0.0
            tZ1 = 0.0
            cs1 = MHDCoordSys.MHDCoordSys(x1, y1, z1,tX1,tY1,tZ1,csSys1)
            r1 = [cs1.x-rmCS.x, cs1.y-rmCS.y, cs1.z - rmCS.z]
            r1 = numpy.array(r1)
            print('r:%s'%r1)
            realDist = (r1[0]**2 + r1[1]**2 + r1[2]**2)**(0.5)
            print('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r1_mag = invRot.apply(r1)
            print('r1_mag:%s'%r1_mag)

            b1Real_mag = rm.getB(r1_mag)
            print('b1Real_mag:%s'%b1Real_mag)
            b1Real_mgtr1 = rot.apply(b1Real_mag)
            print('b1Real_mgtr1:%s'%b1Real_mgtr1)
            b1Real_mgtr1Corr = invRotCorrective.apply(b1Real_mgtr1)
            print('b1Real_mgtr1Corr:%s' % b1Real_mgtr1Corr)
            rCalc_mgtr1 = rm.getR(b1Real_mgtr1Corr)
            print('r1Calc_mgtr1:%s'%rCalc_mgtr1)
            print('r1Calc_magnitude:%s'%(rCalc_mgtr1[0]**2 + rCalc_mgtr1[1]**2 + rCalc_mgtr1[2]**2)**(0.5))

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
            print('r:%s'%r2)
            realDist = (r2[0]**2 + r2[1]**2 + r2[2]**2)**(0.5)
            print('realDistance:%s'%realDist)

            rot = R.from_euler('xyz',[rmCS.thetaX,rmCS.thetaY,rmCS.thetaZ],degrees=True)
            invRot = rot.inv()
            r2_mag = invRot.apply(r2)
            print('r2_mag:%s'%r2_mag)

            b2Real_mag = rm.getB(r2_mag)
            print('b2Real_mag:%s'%b2Real_mag)
            b2Real_mgtr2 = rot.apply(b2Real_mag)
            print('b2Real_mgtr2:%s'%b2Real_mgtr2)
            b2Real_mgtr2Corr = invRotCorrective.apply(b2Real_mgtr2)
            print('b2Real_mgtr2Corr:%s' % b2Real_mgtr2Corr)
            r2Calc_mgtr2 = rm.getR(b2Real_mgtr2Corr)
            print('r2Calc_mgtr1:%s'%r2Calc_mgtr2)
            print('r2Calc_magnitude:%s'%(r2Calc_mgtr2[0]**2 + r2Calc_mgtr2[1]**2 + r2Calc_mgtr2[2]**2)**(0.5))

            #eval positions
            R2R1Real_sys = [cs2.x-cs1.x,cs2.y-cs1.y,cs2.z-cs1.z]
            R2R1RealMag=(R2R1Real_sys[0] ** 2 + R2R1Real_sys[1] ** 2 + R2R1Real_sys[2] ** 2) ** (0.5)
            print('R2R1RealMag:%s' %R2R1RealMag)
            #R2R1Real_sys = rotCorr.apply(R2R1Real_mag)
            print('R2R1Real_sys:%s' % R2R1Real_sys)

            R2R1calc_mag = [r2Calc_mgtr2[0]-rCalc_mgtr1[0],r2Calc_mgtr2[1]-rCalc_mgtr1[1],r2Calc_mgtr2[2]-rCalc_mgtr1[2]]
            R2R1CalcMag = (R2R1calc_mag[0] ** 2 + R2R1calc_mag[1] ** 2 + R2R1calc_mag[2] ** 2) ** (0.5)
            print('R2R1CalcMag:%s' % R2R1CalcMag)
            #might try converting back to global/system coords
            R2R1calc_sys = rotCorrective.apply(R2R1calc_mag)
            print('R2R1Calc_sys:%s' % R2R1calc_sys)


            eVec = [R2R1Real_sys[0]-R2R1calc_sys[0],R2R1Real_sys[1]-R2R1calc_sys[1],R2R1Real_sys[2]-R2R1calc_sys[2]]

            error = R2R1RealMag-R2R1CalcMag
            print('error:%s'%error)
            goalError=1.438244776413855e-05*2.0
            errorFactor = error/goalError
            print('errorFactor:%s'%errorFactor)
            inputs.append(input)
            outputs.append(error)
            out1.append(eVec[0])
            out2.append(eVec[1])
            out3.append(eVec[2])

        for i in range(0,len(inputs)):
            print('%s : %s'%(inputs[i],outputs[i]))

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






tests = [4]

if 2 in tests:
    sys = MHDSystemA(0,0)
    sys.test2()

if 3 in tests:
    sys = MHDSystemA(0,0)
    sys.test3()

if 4 in tests:
    sys = MHDSystemA(0,0)
    sys.test4()