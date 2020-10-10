#!/home/kjeong23/softwares/bin/python3.4
#For a time autocorrelation function, integrate and get exp fit.
#algorithm : integration can be done with Simpson's rule.
#to do exp-fit, first get log of the y value.
#output : screen output with integral, expfit formula.
#file output with log10 of y, fitting line of y.
#input : ACF curve information (t, C(t))

import math
import sys
import numpy
from scipy.integrate import simps
from scipy.optimize import curve_fit
from scipy import stats

def function(x,a,b,sigma): #give x variables and coefficients: calc c0 + c1 cos x. to fit c0,c1
  #c0,c1=coeffs[0],coeffs[1]
  y=a*numpy.exp(-(x-b)*(x-b)/(2*sigma*sigma))
  return y

#def residuals(coeffs,targety,x):
#  return targety-function(x,coeffs)

def main():
  infile=sys.argv[1]
  outfile=open(sys.argv[2],'w')

  octdata=numpy.loadtxt(infile)
  x,targety=octdata[:,0],octdata[:,1]

  #section 1: integral
  intetau=simps(targety,x)
  print("Integrated tau {:11.6f}".format(intetau))

  #section 2: exp fit (linear fit of ln y)
  #before trying linear fit, should first 
  raw_lny=numpy.log(numpy.absolute(targety))
  #raw_log10y=numpy.log10(targety)
  #first, get regional average (5 points)
  #find convergence point of increment : when slope stops changing
  #so |2nd derivative| is close to zero( < 1/100 of max(2nd deriv)), beyond 2nd_deriv max location
  rave_lny=numpy.zeros(len(raw_lny)-2)
  for i in range(len(rave_lny)):
    rave_lny[i]=numpy.average(raw_lny[i:i+3])
  dif_lny,dif2_lny=numpy.zeros(len(rave_lny)-1),numpy.zeros(len(rave_lny)-2)
  for i in range(len(dif_lny)):
    dif_lny[i]=rave_lny[i+1]-rave_lny[i]
  for i in range(len(dif2_lny)):
    dif2_lny[i]=numpy.absolute(dif_lny[i+1]-dif_lny[i])
  #print(dif2_lny)
  j=int(len(dif2_lny)/2)
  maxloc_dif2,max_dif2=numpy.argmax(dif2_lny[:j]),numpy.amax(dif2_lny[:j])
  for i in range(maxloc_dif2,len(dif2_lny)):
    if dif2_lny[i]<max_dif2/1000.0:
      linpoint=i #beginning point of linearity
      break

  rbound=int(len(x)*0.8)
  partx,partlny=x[linpoint:rbound],raw_lny[linpoint:rbound]
  #new_coeffs,cov=curve_fit(function,x,targety,method='lm')  
  #new_coeffs,flag=leastsq(residuals,init_coeffs,args=(waveform_1,x))
  #y=function(x,new_coeffs[0],new_coeffs[1],new_coeffs[2])
  #print("a,b,sigma : ",new_coeffs)
  #print("covariance : ",cov)
  slope, intercept, r_value, p_value, std_err = stats.linregress(partx,partlny)
  print("lnfit slope,intercept,r_value,p_value,std_err,interval_start : \n")
  print("{:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f} \n".format(slope,intercept,r_value,p_value,std_err,x[linpoint]))
  print("lnfit relaxation time: {:11.6f} \n".format(-1.0/slope))

  fitline_original=numpy.exp(slope*x+intercept)
  #for i in range(len(x)):
  #  if i<len(dif2_lny):
  #    outfile.write("{:11.6f} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n".format(x[i],targety[i],fitline_original[i],dif_lny[i],dif2_lny[i]))
  #  else:
  #    outfile.write("{:11.6f} {:11.6f} {:11.6f} 0.0 0.0 \n".format(x[i],targety[i],fitline_original[i]))
  for i in range(len(x)):
    outfile.write("{:11.6f} {:11.6f} {:11.6f} \n".format(x[i],targety[i],fitline_original[i]))

if __name__ == "__main__": main()
