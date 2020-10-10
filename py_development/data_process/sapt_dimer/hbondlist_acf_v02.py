#!/home/kjeong23/softwares/bin/python3.4
# program to calculate hydrogen bonding lifetime.
# algorithm : read H-bond topology index list file.
# maybe not able to construct full h-matrix due to memory issue.
# instead, set a topology index range to evaluate time ACF.
# (h=1 for H-bond, h=0 for non-Hbond for a specific atom triplet)
# By comparing the list of hydrogen bonds of each snapshot, calculate h(t)
# C(t) = (<h(0)h(t)> - <h>^2) / (<h2>-<h>2)
# citation : J. Chem. Phys. 148, 193843 (2018)

import math
import sys
import numpy
import timeit
import mdtraj
import copy

#tstep=1.000 #used ps time unit. 1 snapshot : 1ps.

def autocorr_simple(x,tlimit): # <x(0)x(t)> without considering <x>,<x2>
  #tlimit=int(len(x)/2)
  #a,s2=numpy.mean(x),numpy.var(x)
  zerot,zerotcount=numpy.zeros(tlimit),numpy.zeros(tlimit)
  #zerot[0],zerotcount[0]=1.000,1.000
  if numpy.sum(x)==0: #empty row
    for j in range(0,tlimit):
      zerotcount[j]=len(x)-j
  else:
    for j in range(0,tlimit):
      zerotcount[j]=len(x)-j
      for i in range(len(x)-j):
        zerot[j]+=x[i]*x[i+j]
        #zerotcount[j]+=1
  #ct=numpy.divide(ct,ctcount)

  return zerot,zerotcount

#main fxn
def main():
  #Part to load coordinate file
  hlistfile = open(sys.argv[1],'r')
  outfile = open(sys.argv[2],'w')

  tstep=float(input("How much is the timstep between snapshots in dcd file (in ps)? ex) 10\n"))
  toprange=input("D-H--A triplet index range for multi-cpu calc? ex) 0 10000 (if unnecessary put -1 -1) \n")
  taulimit=int(input("How many snapshot interval to set as upper bound of correlation tau? ex) 500\n"))
  taulimit=taulimit+1
  if toprange=="-1 -1":
    topstart,topend=-1,-1
  else:
    toprsplit=toprange.split()
    topstart,topend=int(toprsplit[0]),int(toprsplit[1])

  start_time=timeit.default_timer()

  #input 1 : load hbond topology index list file
  hblist=[]
  for line in hlistfile:
    hsplit=line.split()
    hsplit=list(map(int,hsplit))
    hblist.append(hsplit)
  nstep=len(hblist)

  elapsed=timeit.default_timer() - start_time
  print('finished hlistfile loading : #snapshots {} time {}'.format(nstep,elapsed))

  #h0ht calculation,collection
  h0ht,h0htcount=numpy.zeros(taulimit,dtype=float),numpy.zeros(taulimit)
  hsum=0
  for i in range(topstart,topend):
    hrow=numpy.zeros(nstep,dtype=int)
    for t in range(nstep):
      if i in hblist[t]:
        hrow[t]=1
    h0ht_one,h0htcount_one=autocorr_simple(hrow,taulimit)
    h0ht+=h0ht_one
    h0htcount+=h0htcount_one
    hsum+=numpy.sum(hrow)
    if i%50==0:
      elapsed=timeit.default_timer() - start_time
      print('hrow correlation complete : index {} time {}'.format(i,elapsed))

  for i in range(len(h0ht)):
    outfile.write('{:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n'.format(i*tstep,h0ht[i],h0htcount[i],hsum,nstep*(topend-topstart)))

  elapsed=timeit.default_timer()-start_time
  print('finished job {}'.format(elapsed))
  outfile.close()
 
if __name__ == "__main__": main()

