#!/usr/bin/python
import random
import os
sourcedir="../program/"
execute="XXZ"
Dim = 3
Latticename = "Pyrochlore"
homedir=os.getcwd()
filelist=os.listdir(sourcedir)
sourcename=[elem for elem in filelist if elem[-3:]=="f90"]
sourcename.sort()
sourcename=sourcename[-1]
os.system("gfortran "+sourcedir+"/"+sourcename+" -O3 -o "+homedir+"/"+execute)
infilepath=homedir+"/infile"
outfilepath=homedir+"/outfile"
inlist=open(homedir+"/inlist","r")
i=0
nblck=int(inlist.readline().split(":")[1])
nsample=int(inlist.readline().split(":")[1])
nsweep=int(inlist.readline().split(":")[1])
ntoss=int(inlist.readline().split(":")[1])
isload=int(inlist.readline().split(":")[1])
seed=-int(random.random()*1000)
if(os.path.exists(infilepath)!=True):
    os.system("mkdir "+infilepath)
if(os.path.exists(outfilepath)!=True):
    os.system("mkdir "+outfilepath)
for eachline in inlist:
        i+=1
        para=eachline.split()
        for j in range(1,int(para[-1])):
            seed-=1
            infile="_in"+str(i)+"_"+str(j)
            outfile="_out"+str(i)+"_"+str(j)
            jobfile="_job"+str(i)+"_"+str(j)+".sh"
            f=open(infilepath+"/"+infile,"w")
            item=para[0:-1]
            item.append(str(ntoss))
            item.append(str(nsample))
            item.append(str(nsweep))
            item.append(str(seed))
            item.append(str(nblck))
            item.append(str(isload))
            item.append("_".join(para[0:-1])+"_"+str(j)+".dat")
            stri=" ".join(item)
            f.write(str(Dim)+"\n")
            f.write("Pyrochlore "+stri)
            f.close()
            elem=str(i)
            os.system("./"+execute+" < "+infilepath+"/"+infile+" > "+outfilepath+"/"+outfile)
