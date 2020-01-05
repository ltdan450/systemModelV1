import MHDCoordSys
import numpy
import scipy

class MHDMagnetometer:

    def __init__(self,cs,system):
        self.cs = cs
        self.system = system
        self.universe = system.universe
        self.mag2mgtrPos = [0.0, 0.0, 0.0]

    ######### normal operation methods#########
    def measure(self):
        print('measure start')
        return self.getBMeasurementFromUniverse()

    #########omnipotent methods#########
    def getBMeasurementFromUniverse(self):
        bField = self.system.universe.getBFieldForMgtr(self)
        #print('bField:%s'%bField)
        return bField