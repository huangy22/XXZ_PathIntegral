#!/usr/bin/env python
import numpy as np
import math
from logger import *
PI=np.pi

class Lattice:
    def __init__(self, Name, Map):
        self.__Map=Map
        self.L=np.array(Map.L)
        self.Name=Name
        try:
            getattr(self, "_Lattice__"+Name)()
        except:
            Assert(False, "Lattice {0} has not been implemented!".format(self.Name))

    def __Pyrochlore(self):
        self.Dim=3
        self.NSublat=4
        self.LatVec=np.array([[0.0,0.5,0.5],
                              [0.5,0.0,0.5],
                              [0.5,0.5,0.0]])
        self.SubLatVec=np.array([[0.0,0.0,0.0],
                                 [0.0,0.25,0.25],
                                 [0.25,0.0,0.25],
                                 [0.25,0.25,0]])

    def GetTetrahedra(self, OnePoint, Sublat, Type):
        ZeroPoint = np.array(OnePoint)-self.SubLatVec[Sublat]
        x = -ZeroPoint[0]+ZeroPoint[1]+ZeroPoint[2]
        y =  ZeroPoint[0]-ZeroPoint[1]+ZeroPoint[2]
        z =  ZeroPoint[0]+ZeroPoint[1]-ZeroPoint[2]
        Vector=[x, y, z]
        Points=[]
        if Type==0:
            for i in range(self.NSublat):
                Points.append(np.array(ZeroPoint+self.SubLatVec[i]))
        elif Type==1:
            for i in range(self.NSublat):
                Points.append(np.array(ZeroPoint-self.SubLatVec[i]))
        return Vector, Points


if __name__=="__main__":
    import weight
    WeightPara={"NSublat": 4, "L":[8, 8, 8],
            "Beta": 0.5, "MaxTauBin":64}
    Map=weight.IndexMap(**WeightPara)
    l=Lattice("Pyrochlore", Map)

    Vector, Points= l.GetTetrahedra([0,0,0], 0, 0)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0,0.5,0.5], 0, 0)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0.5,0,0.5], 0, 0)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0.5,0.5,0], 0, 0)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0,1,0], 0, 0)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0.5,1,0.5], 0, 0)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0,0.5,0.5], 0, 1)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0.5,0,0.5], 0, 1)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0.5,0.5,0], 0, 1)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0,1,0], 0, 1)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0.5,1,0.5], 0, 1)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([1.0,0.5,0.5], 0, 1)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'

    Vector, Points= l.GetTetrahedra([0.5,0.5,1.0], 0, 1)
    for i in range(l.NSublat):
        print '[', Points[i], ',',  Vector, ',', i, '],'
