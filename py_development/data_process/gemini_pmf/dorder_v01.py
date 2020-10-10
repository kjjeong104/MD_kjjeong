#!/home/kjeong23/softwares/bin/python3.4
# Directional order parameter calculator
# algorithm: define two vectors and interval to calculate ensemble average
# then calculate cosine, finally calc order parameter S (smectic order p)
# get atomindex list of solute acceptor -> read ndxfile(1step)(3column:donorO-H-acceptor) 
# -> starting from solute acceptor, assign all water molecules' Hbond network layers
# -> count #Hbond for all water molecules, and average by layers -> output
# 3 arguments : vector1file vector2file output (also need input for block settings)
# caution : the input file should be .xvg format produced by gromacs g_dist module

import math
import sys
import numpy
import timeit

def main():
  #Load files
  v1file = open(sys.argv[1],'r')
  v2file = open(sys.argv[2],'r')
  outfile = open(sys.argv[3],'w') #output file for order parameter value of the block

  start_time=timeit.default_timer()
  RADDEG=57.29577951
  nskip=22 #number of lines we should skip in xvg file
  #ask for block settings
  nblocks=int(input("Number of time blocks for block averaging? ex) 3\n"))
  initt=float(input("Initial time of traj in ps? ex) 7000\n"))
  blength=float(input("Length of one time block in ps? ex) 15000\n"))

  S=numpy.zeros(nblocks) #smectic order parameter
  nstep=numpy.zeros(nblocks) ## of steps for each block
  nline=0

  #load vector files and calc values on the fly
  while nline<nskip:
    line1=v1file.readline()
    line2=v2file.readline()
    nline+=1

  while True:
    line1=v1file.readline()
    line2=v2file.readline()
    if (""==line1 or ""==line2):
      break;
    set1=line1.split()
    set2=line2.split()
    info1=numpy.array([float(x) for x in set1])
    info2=numpy.array([float(x) for x in set2])
    #print(info1,info2)
    if (info1[0]==info2[0] and info1[0]>=initt):
      b=int( (info1[0]-initt)/blength ) #block index
      if b in range(nblocks): #calculate
        cosd=numpy.dot(info1[2:5],info2[2:5])/(info1[1]*info2[1])
        #print(cosd, numpy.arccos(cosd), info1[0], b)
        smalls=(3.0*cosd*cosd-1.0)/2.0 #order parameter of one incident
        nstep[b]+=1.0
        S[b]+=smalls
        #print(info1[2:5],info2[2:5],numpy.dot(info1[2:5],info2[2:5]),cosd,smalls)

  elapsed=timeit.default_timer() - start_time
  print('collection complete time {}'.format(elapsed))

 #printing section
  for i in range(nblocks):
    S[i] = S[i]/nstep[i]
    outfile.write('{:8.4f} '.format(S[i]))

  print(nstep)

  v1file.close()
  v2file.close()
  outfile.close()

if __name__ == "__main__": main()


