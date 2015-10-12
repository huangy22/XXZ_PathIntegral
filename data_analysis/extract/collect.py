#!/usr/bin/python
import os
n=4
Beta = 5.0
L = 4
path = "../../Job6/"
os.system("rm *s_sqa*.dat")
os.system("rm corr_*.txt")
os.system("rm static_cor*.txt")
os.system("gfortran ../../program/extdata.f90 -o ./extdata")
os.system("cp "+path+"mid_hs_sqa0_*.txt .")
os.system("gfortran ../../program/extcorr.f90 -o ./extcorr")
os.system("cp "+path+"corr_k*.txt .")
os.system("cp "+path+"static_corr_*.txt .")

for i in range(0,n):
    order="echo "+str(i)+" "+str(n)+" >temp.dat"
    os.system(order)
    order="./extdata < temp.dat"
    os.system(order)

order="echo "+str(n)+" "+str(L)+" "+str(Beta)+" >tempcorr.dat" 
os.system(order)
order="./extcorr < tempcorr.dat"
os.system(order)

os.system("cp static_corr.txt "+path)
os.system("cp corr_k1.txt "+path)
os.system("cp corr_k2.txt "+path)
os.system("cp hs_sqa0.dat "+path)
os.system("rm temp.dat")
os.system("rm tempcorr.dat")

