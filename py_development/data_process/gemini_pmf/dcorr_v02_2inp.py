#!/home/kjeong23/softwares/bin/python3.4
# Directional autocorrelation function calculator
# algorithm: define vector, get production interval info, store all vectors in time array
# then calculate self-correlating cosine, order parameter S, bin in time corr fxn array
# averaging 
# v02 : can accept 2 input files and average them (tailored for gemini 2body problem), exp fit
# 3 arguments : vectorfile1 vectorfile2 output (also need input for production run settings)
# caution : the input file should be .xvg format produced by gromacs g_dist module

import math
import sys
import numpy
from scipy.optimize import curve_fit
import timeit

def expfunc(x,a,b):
  return a*numpy.exp(-b*x)

def main():
  #Load files
  v1file = open(sys.argv[1],'r')
  v2file = open(sys.argv[2],'r')
  outfile = open(sys.argv[3],'w') #output file for order parameter value of the block

  start_time=timeit.default_timer()
  RADDEG=57.29577951
  nskip=22 #number of lines we should skip in xvg file
  #ask for block settings
  #nblocks=int(input("Number of time blocks for block averaging? ex) 3\n"))
  initt=float(input("Initial time of traj in ps? ex) 7000\n"))
  endt=float(input("Final time of traj in ps? ex) 52000\n"))
  tstep=float(input("time between data in ps? ex) 10\n"))
  corrlim=float(input("Correlation time range limit in ps? putting 0 will autoset ex) 1000\n"))

  totstep=int( (endt-initt)/tstep )
  corstep=int(totstep/2) #number of tau intervals that actually count time corr
  if corrlim!=0:
    corstep=int(corrlim/tstep)+1
  Pet=numpy.zeros(corstep) #end-to-end vector time corr fxn
  vecf1=numpy.zeros((totstep,4))
  vecf2=numpy.zeros((totstep,4))
  #S=numpy.zeros(nblocks) #smectic order parameter
  nstep=numpy.zeros(corstep) ## of steps for each bin in time corr fxn
  nline=0

  #load vector files and store vector information
  while nline<nskip:
    line1=v1file.readline()
    line2=v2file.readline()
    nline+=1
  nline=0
  while True:
    line1=v1file.readline()
    if ""==line1:
      break;
    set1=line1.split()
    info1=numpy.array([float(x) for x in set1])
    #print(info1,info2)
    if (info1[0]>=initt and info1[0]<endt):
      vecf1[nline]=info1[1:5]
      nline+=1
  nline=0
  while True:
    line2=v2file.readline()
    if ""==line2:
      break;
    set2=line2.split()
    info2=numpy.array([float(x) for x in set2])
    #print(info1,info2)
    if (info2[0]>=initt and info2[0]<endt):
      vecf2[nline]=info2[1:5]
      nline+=1

  #corr function calculation (double loop)
  for i in range(totstep): #initial time
    for j in range(corstep): #correlation time
      #print(i,j)
      if j < totstep-i:
        v1,v2=vecf1[i],vecf1[i+j]
        cosd=numpy.dot(v1[1:4],v2[1:4])/(v1[0]*v2[0])
        smalls=(3.0*cosd*cosd-1.0)/2.0 #order parameter at the incident
        nstep[j]+=1.0
        Pet[j]+=smalls
        v1,v2=vecf2[i],vecf2[i+j]
        cosd=numpy.dot(v1[1:4],v2[1:4])/(v1[0]*v2[0])
        smalls=(3.0*cosd*cosd-1.0)/2.0 #order parameter at the incident
        nstep[j]+=1.0
        Pet[j]+=smalls

  #time array
  tarray=numpy.zeros(corstep)
  for j in range(corstep):
    tarray[j]=j*tstep

  #printing section
  for j in range(corstep):
    Pet[j] = Pet[j]/nstep[j]
    outfile.write('{:8.4f} {:8.4f} \n'.format(tarray[j],Pet[j]) )

  print(nstep)

  #exponential fitting
  popt,pcov=curve_fit(expfunc,tarray,Pet)
  relaxt=1/popt[1]
  #print(popt)

  outfile.write('Relaxtime_in_ps: {:8.4f} \n'.format(relaxt) )

  v1file.close()
  v2file.close()
  outfile.close()

if __name__ == "__main__": main()

