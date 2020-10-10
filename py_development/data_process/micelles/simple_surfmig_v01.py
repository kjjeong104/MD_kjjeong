#!/home/kjeong23/softwares/bin/python3.4
# micelle kinetics analysis program
# focus on surfactant migration. track how individual surfactant is migrated.
# input file : merged lattice point assignment summary file for entire time period.(reassigned)

import math
import sys
import numpy
import timeit

mthres=6 #threshold aggregation number of an independent micelle : 6
def comparison(a1,a2): #getting retention fraction of micelle during change a1->a2
  coin=0
  for x in a1:
    if x in a2:
      coin+=1
  retf=coin/len(a1)
  return retf

def infosplit(line):
  memline=line[96:]
  memline=memline.replace(",","")
  memsplit=memline.split()
  memsplit=[int(x) for x in memsplit]
  stsplit=line[6:18].split()
  site,time=int(stsplit[0]),float(stsplit[1])
  infoline=[site]+[time]+memsplit
  return infoline

#main fxn
def main():
  #Load input files
  morfile = open(sys.argv[1],'r')
  outfile = open(sys.argv[2],'w') #output file for kinetics summary
  #matrixfile = open(sys.argv[3],'w')

  #set parameters if needed
  tstep=float(input("Time difference between snapshots (in ns)? ex) 0.02 (stop if irregular)\n"))
  nnmic=int(input("Most frequent number of micelles desired? ex) 30 or 8\n"))
  nsurf=int(input("Number of surfactants? ex) 973\n"))

  start_time=timeit.default_timer()
  #read the clusterization result file, store all data in arrays
  #core array:3d list - [step][micelle][member]
  i,nstep=0,0
  biglist,steplist,nmiclist=[],[],[]
  for line in morfile:
    if line[0:4]=="step":#step initiating line
      ltsplit=line.split()
      nmcl1=int(ltsplit[3]) #number of mcls, temporary use
      nmiclist.append(nmcl1)
      i=0
      nstep+=1
    else: #collect
      i+=1 
      infolist=infosplit(line)
      steplist.append(infolist)
      if i==nmcl1: #1 timestep ended
        biglist.append(steplist)
        steplist=[]
  elapsed=timeit.default_timer() - start_time
  print('finished reading cluster file, {} snapshots detected, time {:11.4f}'.format(nstep,elapsed))
  print('data size {} bytes'.format(sys.getsizeof(biglist)))

  #construct surfactant location matrix
  surfloc=numpy.zeros((nsurf,nstep))
  i=0
  for step1 in biglist:
    for row in step1:
      mem=row[2:]
      for x in mem:
        surfloc[int(x)][i]=int(row[0])
    i+=1

  totmig_surf,totacc_mig,totdon_mig=numpy.zeros(nsurf),numpy.zeros(nnmic),numpy.zeros(nnmic)
  #surfactant migration output printing section
  outfile.write("surf# initmic# destin# initt final\n")
  for i in range(nsurf):
    recent=surfloc[i][0] #residing micelle at step0
    recstep=0
    for j in range(1,nstep):
      newmic=surfloc[i][j]
      if newmic!=0:
        if newmic!=recent: #detected surfactant migration
          # format: surf# initmic# destination# initt finalt
          outfile.write('{:4}  {:3}   {:3} {:8.3f} {:8.3f}\n'.format(i,recent,newmic,recstep*tstep,j*tstep))
          totmig_surf[i]+=1
          totacc_mig[int(newmic)-1]+=1
          totdon_mig[int(recent)-1]+=1
        recent=surfloc[i][j]
        recstep=j
  outfile.write("==========end of events========below:total counts\n")
  for i in range(nsurf):
    outfile.write('surf# {:4} migration {:5}\n'.format(i,totmig_surf[i]))
  for i in range(nnmic):
    outfile.write('mic# {:3} acc {:5} don {:5}\n'.format(i+1,totacc_mig[i],totdon_mig[i]))

  morfile.close()
  outfile.close()

if __name__ == "__main__": main()

