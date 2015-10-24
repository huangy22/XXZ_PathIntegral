#!/usr/bin/python
import os
n=200
Beta = 500
L = 4
Jx=0.09
path = "../"
os.system("rm *s_sqa*.dat")

os.system("gfortran ./extdata.f90 -o ./extdata")
os.system("gfortran ./extcorr.f90 -o ./extcorr")
os.system("gfortran ./extract.f90 -o ./extract")
os.system("cp "+path+"mid_hs_sqa0_*.txt .")

for i in range(0,n):
    order="echo "+str(i)+" "+str(n)+" >temp.dat"
    os.system(order)
    order="./extdata < temp.dat"
    os.system(order)

order="./extract"
os.system(order)

os.system("rm mid_hs_sqa*.dat")
os.system("cp hs_sqa0.dat "+path)
os.system("cp ext_hs_sqa0.dat "+path)
