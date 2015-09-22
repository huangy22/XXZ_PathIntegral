#!/usr/bin/python
import os
n=2
temp="temp.dat"
os.system("rm hs_sqa0.dat")
os.system("gfortran ../../program/extdata.f90 -o ./extdata")
os.system("cp ../../Job1/mid_hs_sqa0_*.txt .")
for i in range(0,n):
    order="echo "+str(i)+" "+str(n)+" >temp.dat"
    os.system(order)
    order="./extdata < temp.dat"
    os.system(order)

os.system("rm temp.dat")

