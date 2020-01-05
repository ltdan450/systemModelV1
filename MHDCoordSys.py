import numpy
import scipy

class MHDCoordSys:
    def __init__(self, xP=0, yP=0, zP=0, thetaXP=0, thetaYP=0, thetaZP=0, refCS=0):
        self.x = xP
        self.y = yP
        self.z = zP
        self.thetaX = thetaXP
        self.thetaY = thetaYP
        self.thetaZ = thetaZP
        self.rot = 'xyz'
        self.refCS = refCS

        self.childrenCS = []

        if refCS:
            self.refCS.addChildCS(self)

    def transformVectFromCS(self, vect, fromCS):
        print('TBD')

    def addChildCS(self,childCS):
        if childCS not in self.childrenCS:
            self.childrenCS.append(childCS)