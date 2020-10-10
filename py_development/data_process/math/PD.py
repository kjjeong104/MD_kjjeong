#!/home/kjeong23/softwares/bin/python3.4
#probability distribution function
#two inputs : inputfile name, outputfile name

import math
import sys
import numpy

def main():
  infile=sys.argv[1]
  outfile=open(sys.argv[2],'w')

  brange = input("Min,Max and Grid size for Prob.dist? ex) 0 80 1\n")
  integ = float(input("Integral of prob.dist after normalization? ex) 1.0 \n"))
  bset=brange.split()
  bmin=float(bset[0])
  bmax=float(bset[1])
  bbin=float(bset[2])
  nbin=int((bmax-bmin)/bbin)+1 #include right endpoint

  #b=numpy.empty(nbin)
  x=numpy.loadtxt(infile) #read raw data. len(x) will serve as total count
  count=len(x)
  pd=numpy.zeros(nbin)
  #for i in range(nbin):
  #  b[i]=bmin+bbin*i

  for dat in x:
    index=int((dat-bmin+bbin/2.0)/bbin)
    if index<nbin:
      pd[index]+=1   

  a,s2=numpy.mean(x),numpy.var(x)
  print("ave,var is ",a,s2)

  #normalize
  pd=pd/(count*bbin) * integ
  for i in range(nbin):
    outfile.write('{:15.8f} {:15.8f}\n'.format(b[i],pd[i]))
  outfile.close()

if __name__ == "__main__": main()

