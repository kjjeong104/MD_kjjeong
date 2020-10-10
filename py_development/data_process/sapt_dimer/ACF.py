#!/usr/bin/env python
#autocorrelation function
#two inputs : inputfile name, outputfile name

import math
import sys
import numpy

def main():
  infile=sys.argv[1]
  outfile=open(sys.argv[2],'w')

  x=numpy.loadtxt(infile)

  tlimit=int(len(x)/2)
  a,s2=numpy.mean(x),numpy.var(x)
  ct,ctcount=numpy.zeros(tlimit),numpy.zeros(tlimit)
  ct[0],ctcount[0]=1.000,1.000
  #print(tlimit,a,s2)

  for j in range(1,tlimit):
    for i in range(len(x)-j):
      ct[j]+=(x[i]-a)*(x[i+j]-a)/s2
      ctcount[j]+=1

  ct=numpy.divide(ct,ctcount)
  #print(ctcount) 
  for x in ct:
    outfile.write('{:11.6f}\n'.format(x))
  outfile.close()

if __name__ == "__main__": main()

