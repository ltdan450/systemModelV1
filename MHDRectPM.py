import numpy
import math
import random
import MHDCoordSys
from scipy.spatial.transform import Rotation as R

verbose = 0
verbose = 0
def log(input):
    if verbose:
        print(input)

class MHDRectPM:
    def __init__(self,csIn,LxIn,LyIn,LzIn,JIn):
        #self.Lx = LxIn
        #self.Ly = LyIn
        #self.Lz = LzIn
        self.J = JIn
        if 1 == 1:
            self.a = LxIn
            self.b = LyIn
            self.h = LzIn
        else:
            self.a = LyIn
            self.b = LzIn
            self.h = LxIn
        self.muZero = 4.0 * 3.141592653 * 1.0e-7
        self.K = self.muZero * self.J / (4.0 * 3.1415926535)

        self.x = -.14
        self.y = 0.001
        self.z = -0.001

        self.dBxDx = .1
        self.dByDx = .1
        self.dBxDy = .1
        self.dByDy = .1

        self.sigmaBx = 0.293
        self.sigmaBy = 0.164
        self.sigmaBz = 0.207

        xOld = .1
        yOld = .1

        self.xMin = -1e5
        self.xMax = 1e5
        self.yMin = -1e5
        self.yMax = 1e5
        self.zMin = 0
        self.zMax = 1e5

        self.cs = csIn

        self.orientationType = 1

        self.verbose = 0

        self.posVect = numpy.array([[1.0],[2.0],[1.5]],numpy.double)

    def gammaP(self,g1,g2,g3,zP0):
        gOut = 1.0

        numerator = math.sqrt(g1 * g1 + g2 * g2 + (g3 - zP0) * (g3 - zP0)) - g2
        denominator = math.sqrt(g1 * g1 + g2 * g2 + (g3 - zP0) * (g3 - zP0)) + g2
        gOut = math.log(numerator / denominator)

        return gOut

    def phiP(self,p1,p2,p3,zP0):
        numerator = p1 * (p3 - zP0)
        denominator = p2 * math.sqrt(p1 * p1 + p2 * p2 + (p3 - zP0) * (p3 - zP0))
        if (denominator != 0.0):
            pOut = math.atan(numerator / denominator)
        else:
            pOut = 0.0
        return pOut

    def BxP(self,xP,yP,zP):
        BxPH = self.gammaP(self.a - xP, yP, zP, self.h) + self.gammaP(self.a - xP, self.b - yP, zP, self.h) - self.gammaP(xP, yP, zP, self.h) - self.gammaP(xP,self.b - yP,zP,self.h)
        BxP0 = self.gammaP(self.a - xP, yP, zP, 0.0) + self.gammaP(self.a - xP, self.b - yP, zP, 0.0) - self.gammaP(xP, yP, zP, 0.0) - self.gammaP(xP, self.b - yP, zP, 0.0)
        BxPout = - self.K / 2.0 * (BxPH - BxP0)
        return BxPout

    def ByP(self,xP,yP,zP):
        ByPH = self.gammaP(self.b - yP,xP,zP,self.h)
        ByPH+= self.gammaP(self.b - yP,self.a - xP,zP,self.h)
        ByPH-= self.gammaP(yP,xP,zP,self.h)
        ByPH-= self.gammaP(yP,self.a - xP,zP,self.h)

        ByP0 = self.gammaP(self.b - yP,xP,zP,0.0)
        ByP0+= self.gammaP(self.b - yP,self.a - xP,zP,0.0)
        ByP0-= self.gammaP(yP,xP,zP,0.0)
        ByP0-= self.gammaP(yP,self.a - xP,zP,0.0)

        ByPout = - self.K / 2.0 * (ByPH - ByP0)
        return ByPout

    def BzP(self,xP,yP,zP):
        BzPH = self.phiP(yP,self.a - xP,zP,self.h)          
        BzPH+=self.phiP(self.b - yP,self.a - xP,zP,self.h)
        BzPH+=self.phiP(xP,self.b - yP,zP,self.h)    
        BzPH+=self.phiP(self.a - xP,self.b - yP,zP,self.h)     
        BzPH+=self.phiP(self.b - yP,xP,zP,self.h)
        BzPH+=self.phiP(yP,xP,zP,self.h)       
        BzPH+=self.phiP(self.a - xP,yP,zP,self.h)        
        BzPH+=self.phiP(xP,yP,zP,self.h)

        BzP0 = self.phiP(yP,self.a - xP,zP,0.0)
        BzP0+=self.phiP(self.b - yP,self.a - xP,zP,0.0)
        BzP0+=self.phiP(xP,self.b - yP,zP,0.0)
        BzP0+=self.phiP(self.a - xP,self.b - yP,zP,0.0)
        BzP0+=self.phiP(self.b - yP,xP,zP,0.0)
        BzP0+=self.phiP(yP,xP,zP,0.0)
        BzP0+=self.phiP(self.a - xP,yP,zP,0.0)
        BzP0+=self.phiP(xP,yP,zP,0.0)

        BzPout = - self.K * (BzPH - BzP0)
        return BzPout

    def getUpdatedJacobian(self,xP,yP,zP,delta):
        # Xprime B field partial derivatives
        B0 = self.BxP(xP,yP,zP)
        B1 = self.BxP(xP+delta, yP, zP)
        P_BxP_P_xP = (B1-B0)/delta

        B0 = self.BxP(xP, yP, zP)
        B1 = self.BxP(xP, yP+delta, zP)
        P_BxP_P_yP = (B1 - B0) / delta

        B0 = self.BxP(xP, yP, zP)
        B1 = self.BxP(xP, yP, zP+delta)
        P_BxP_P_zP = (B1 - B0) / delta

        # Yprime B field partial derivatives
        B0 = self.ByP(xP, yP, zP)
        B1 = self.ByP(xP+delta, yP, zP)
        P_ByP_P_xP = (B1 - B0) / delta

        B0 = self.ByP(xP, yP, zP)
        B1 = self.ByP(xP, yP+delta, zP)
        P_ByP_P_yP = (B1 - B0) / delta

        B0 = self.ByP(xP, yP, zP)
        B1 = self.ByP(xP, yP, zP+delta)
        P_ByP_P_zP = (B1 - B0) / delta

        # Zprime B field partial derivatives
        B0 = self.BzP(xP, yP, zP)
        B1 = self.BzP(xP+delta, yP, zP)
        P_BzP_P_xP = (B1 - B0) / delta

        B0 = self.BzP(xP, yP, zP)
        B1 = self.BzP(xP, yP+delta, zP)
        P_BzP_P_yP = (B1 - B0) / delta

        B0 = self.BzP(xP, yP, zP)
        B1 = self.BzP(xP, yP, zP+delta)
        P_BzP_P_zP = (B1 - B0) / delta

        jac = numpy.array([[P_BxP_P_xP, P_BxP_P_yP, P_BxP_P_zP],[P_ByP_P_xP, P_ByP_P_yP, P_ByP_P_zP],[P_BzP_P_xP, P_BzP_P_yP, P_BzP_P_zP]])

        return jac

    def updatePosition(self,Bx,By,Bz,guessX,guessY,guessZ):
        guessXP = guessX#*-1.0
        guessYP = guessY#*-1.0
        guessZP = guessZ
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
        errorThreshold = 0.000001

        log('begin nonlinear equation solve')

        while (trialPos<=10000):
            log('trialPos loop begin:%i'%trialPos)
            iterations = 0
            while iterations<10:
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
                log('error norm:%f'%BPErrorNorm)

                x = posPrime[0, 0]
                y = posPrime[1, 0]
                z = posPrime[2, 0]

                if (BPErrorNorm < errorThreshold):

                    validPosFlag = 1
                    break                                       #break iterations loop
                iterations+=1
            #proceed
            log('validPosFlag:%i'%validPosFlag)
            if validPosFlag == 1:
                guessXP = posPrime[0,0]
                guessYP = posPrime[1,0]
                guessZP = posPrime[2,0]

                self.x = guessXP
                self.y = guessYP
                self.z = guessZP

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

    def getBP(self,xyzIn):
        bPMag = numpy.array([self.BxP(xyzIn[0],xyzIn[1],xyzIn[2]),self.ByP(xyzIn[0],xyzIn[1],xyzIn[2]),self.BzP(xyzIn[0],xyzIn[1],xyzIn[2])], numpy.double)
        return bPMag

#Test
if 1 == 0:
    #standard home depot magnet used
    #zP = -x
    #xP = -y
    #yP = z

    if 1 == 0:
        Lx = .047
        Ly = .0215
        Lz = .0097
        testMag = MHDRectPM(0.047,0.0215,0.0097,1.25e11)

        xP = Lx/2.0 + .2
        yP = Ly/2.0+ .2
        zMin = 0.0
        zMax = 0.2

        for i in range(0,5,1):
            z = zMin + (zMax-zMin)/4 * i
            print ('test value BxP:%f ByP:%s BzP:%s z:%s'%(testMag.BxP(xP,yP,z),testMag.ByP(xP,yP,z),testMag.BzP(xP,yP,z),z))

    if 1 == 1:
        testMag = MHDRectPM(0.047, 0.0215, 0.0097, 1.25e11)
        bxMeasured = -38.1
        byMeasured = 33.4
        bzMeasured = -.505309
        bM = numpy.array([bxMeasured, byMeasured, bzMeasured], numpy.double)

        xMeasured = -.1719
        yMeasured = 0.06330
        zMeasured = .0030  # (.0456/10)
        posM = numpy.array([xMeasured, yMeasured, zMeasured], numpy.double)
        print('posM:%s' % posM)

        globalToPrime = ['yz', [90, 90]]
        r = R.from_euler(globalToPrime[0], globalToPrime[1], degrees=True)

        posPM = r.apply(posM)
        print('posPM:%s' % posPM)

        bPM = r.apply(bM)
        print('bM:%s' % bM)
        print('bPM:%s' % bPM)

        bxPMeasured = bPM[0]
        byPMeasured = bPM[1]
        bzPMeasured = bPM[2]

        xPMeasured = posPM[0]
        yPMeasured = posPM[1]
        zPMeasured = posPM[2]

        res = testMag.updatePosition(bxPMeasured,byPMeasured,bzPMeasured,xPMeasured,yPMeasured,zPMeasured)

        print('result:%s'%res)

        print('actual xPM:%s yPM:%s zPM:%s'%(xPMeasured,yPMeasured,zPMeasured))
        print('calc   xPM:%s yPM:%s zPM:%s'%(res[0],res[1],res[2]))

        xP = numpy.array([-.72,.005,.126],numpy.double)

        #bXPCalc = numpy.array([testMag.BxP(xP[0],xP[1],xP[2]),testMag.ByP(xP[0],xP[1],xP[2]),testMag.BzP(xP[0],xP[1],xP[2])], numpy.double)
        #print('bXPCalc:%s'%bXPCalc)


        #print(testMag.BzP(0.0,0.0))








#scratch secion
if 1 == 0:
    arr = numpy.array([[1.0],[2.0],[1.5]],numpy.double)

    print(arr)