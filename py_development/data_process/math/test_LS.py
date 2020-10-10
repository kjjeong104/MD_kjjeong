#!/home/kjeong23/softwares/bin/python3.4
#least square fitting
#two inputs : inputfile name, outputfile name

import math
import sys
import numpy
from scipy.optimize import curve_fit

def function(x,c0,c1): #give x variables and coefficients: calc c0 + c1 cos x. to fit c0,c1
  #c0,c1=coeffs[0],coeffs[1]
  y=c1*numpy.cos(x)+c0
  return y

def residuals(coeffs,targety,x):
  return targety-function(x,coeffs)

def main():
  infile=sys.argv[1]

  odata=numpy.loadtxt(infile) #original data set to fit
  #print(odata)
  x=odata[:,0]
  targety=odata[:,1]

  #init_coeffs=numpy.array([1,1],dtype=float) #true solution should be [-50,150]

  new_coeffs,cov=curve_fit(function,x,targety,method='lm')  
  #new_coeffs,flag=leastsq(residuals,init_coeffs,args=(waveform_1,x))
  print(new_coeffs)

if __name__ == "__main__": main()

