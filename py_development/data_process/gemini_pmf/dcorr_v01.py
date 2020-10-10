#!/home/kjeong23/softwares/bin/python3.4
# Directional autocorrelation function calculator
# algorithm: define vector, get production interval info, store all vectors in time array
# then calculate self-correlating cosine, order parameter S, bin in time corr fxn array
# averaging 
# 2 arguments : vectorfile output (also need input for production run settings)
# caution : the input file should be .xvg format produced by gromacs g_dist module

import math
import sys
import numpy
import timeit

def main():
  #Load files
  v1file = open(sys.argv[1],'r')
  outfile = open(sys.argv[2],'w') #output file for order parameter value of the block

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
  vec=numpy.zeros((totstep,4))
  #S=numpy.zeros(nblocks) #smectic order parameter
  nstep=numpy.zeros(corstep) ## of steps for each bin in time corr fxn
  nline=0

  #load vector files and store vector information
  while nline<nskip:
    line1=v1file.readline()
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
      vec[nline]=info1[1:5]
      nline+=1

  #corr function calculation (double loop)
  for i in range(totstep): #initial time
    for j in range(corstep): #correlation time
      #print(i,j)
      if j < totstep-i:
        v1,v2=vec[i],vec[i+j]
        cosd=numpy.dot(v1[1:4],v2[1:4])/(v1[0]*v2[0])
        smalls=(3.0*cosd*cosd-1.0)/2.0 #order parameter at the incident
        nstep[j]+=1.0
        Pet[j]+=smalls

  #printing section
  for j in range(corstep):
    Pet[j] = Pet[j]/nstep[j]
    outfile.write('{:8.4f} {:8.4f} \n'.format(j*tstep,Pet[j]) )

  print(nstep)

  v1file.close()
  outfile.close()

if __name__ == "__main__": main()


