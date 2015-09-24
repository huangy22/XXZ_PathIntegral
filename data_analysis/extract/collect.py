#!/usr/bin/python
import os
n=50
Beta = 5.0
L = 4
os.system("rm hs_sqa0.dat")
os.system("rm corr_frequency.txt")
os.system("rm static_corr.txt")
os.system("gfortran ../../../program/extdata.f90 -o ./extdata")
os.system("cp ../../mid_hs_sqa0_*.txt .")
os.system("gfortran ../../../program/extcorr.f90 -o ./extcorr")
os.system("cp ../../corr_frequency_*.txt .")
os.system("cp ../../static_corr_*.txt .")

for i in range(0,n):
    order="echo "+str(i)+" "+str(n)+" >temp.dat"
    os.system(order)
    order="./extdata < temp.dat"
    os.system(order)

order="echo "+str(n)+" "+str(L)+" "+str(Beta)+" >tempcorr.dat" 
os.system(order)
order="./extcorr < tempcorr.dat"
os.system(order)

os.system("rm temp.dat")
os.system("rm tempcorr.dat")

