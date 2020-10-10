#!/home/kjeong23/softwares/bin/python3.4
#fitting to gaussian curve
#two inputs : inputfile name, outputfile name

import math
import sys
import numpy
from scipy.optimize import curve_fit

def function(x,a,b,sigma): #give x variables and coefficients: calc c0 + c1 cos x. to fit c0,c1
  #c0,c1=coeffs[0],coeffs[1]
  y=a*numpy.exp(-(x-b)*(x-b)/(2*sigma*sigma))
  return y

#def residuals(coeffs,targety,x):
#  return targety-function(x,coeffs)

def main():
  infile=sys.argv[1]
  outfile=open(sys.argv[2],'w')

  odata=numpy.loadtxt(infile) #original data set to fit
  #print(odata)
  x=odata[:,0]
  targety=odata[:,1]

  #init_coeffs=numpy.array([1,1],dtype=float) #true solution should be [-50,150]

  new_coeffs,cov=curve_fit(function,x,targety,method='lm')  
  #new_coeffs,flag=leastsq(residuals,init_coeffs,args=(waveform_1,x))
  y=function(x,new_coeffs[0],new_coeffs[1],new_coeffs[2])
  print("a,b,sigma : ",new_coeffs)
  print("covariance : ",cov)

  for i in range(len(x)):
    outfile.write("{:11.6f} {:11.6f} {:11.6f}\n".format(x[i],y[i],targety[i]))

if __name__ == "__main__": main()

