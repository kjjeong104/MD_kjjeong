#!/home/kjeong23/softwares/bin/python3.4
# Simple text parsing based program to convert
# voronoi inclusion data file into voronoi inclusion density data file
# 3 arguments : incfile volfile output

import math
import sys
import numpy
import timeit

def main():
  #Load files
  incfile = open(sys.argv[1],'r')
  volfile = open(sys.argv[2],'r')
  outfile = open(sys.argv[3],'w') #output file for revised voronoi analysis table

  start_time=timeit.default_timer()
  sindex=-1
  #load incfile, volfile simultaneously. Can read line by line
  for line1 in incfile:
    if line1[0:3]=='mic': #label line
      sindex+=1
      if (sindex%50)==0:
        elapsed=timeit.default_timer() - start_time
        print('step {} started time {}'.format(sindex,elapsed))
      split=line1.split()
      outfile.write('mic# agN    Rg       Asp     sv2a     sv2b     sv2c  \
   waterd    TMAd     massd   chgd   step {} numofmcls {}\n'.format(split[12],split[14])) 
      line2=volfile.readline() #skip volfile line
    else: #data line
      split1=line1.split()
      line2=volfile.readline()
      split2=line2.split()
      counts=numpy.array([float(x) for x in split1])
      for i in range(7,11):
        counts[i]/=float(split2[1]) #dividing 'counts' by voronoi polyhedron volume
      outfile.write('{:3} {:3} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n'\
.format(int(split1[0]),int(split1[1]),counts[2],counts[3],counts[4],counts[5],counts[6],counts[7],counts[8],counts[9],counts[10]))

  incfile.close()
  volfile.close()
  outfile.close()

if __name__ == "__main__": main()
