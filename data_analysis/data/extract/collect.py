#!/usr/bin/python
import os
n=400
Beta = 500
L = 4
Jx=0.09
path = "../"
os.system("rm *s_sqa*.dat")

os.system("gfortran ./extdata.f90 -o ./extdata")
os.system("gfortran ./extcorr.f90 -o ./extcorr")
os.system("gfortran ./extract.f90 -o ./extract")
os.system("cp "+path+"mid_hs_sqa0_*.txt .")
os.system("cp "+path+"corr_k*.txt .")
os.system("cp "+path+"static_corr*.txt .")

for i in range(0,n):
    order="echo "+str(i)+" 64 >temp.dat"
    os.system(order)
    order="./extdata < temp.dat"
    os.system(order)

order="./extract"
os.system(order)

order="echo "+str(n)+" "+str(L)+" "+str(Beta)+" >tempcorr.dat"
os.system(order)
order="./extcorr < tempcorr.dat"
os.system(order)


os.system("cp hs_sqa0.dat "+path)
os.system("cp ext_hs_sqa0.dat "+path)
os.system("cp static_corr.txt "+path)
os.system("cp corr_k*.txt "+path)

os.system("rm temp*.dat")
os.system("rm mid_*.txt")
os.system("rm static_corr_1_*.txt")
os.system("rm corr_k*_1_*.txt")
