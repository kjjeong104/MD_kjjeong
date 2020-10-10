#!/home/kjeong23/softwares/bin/python3.4

import math
import numpy
# Newton-Raphson method practice

#Variables section
thres=0.00001  #threshold
target=0      #target value for the function
guess=0     #initial guess clue for x
dx=0.05    #increment of x to measure slope of f(x). search for both direction
damp=1  #damping factor

#solve equation f(x)=c, when c is the target.
# for this case, let f(x)=x^3-2x+2

def f(x):
  value=0.002*x**2-0.007*x-0.0005
  return value

x=guess
y=f(x)
iter=0

while abs(y-target) > thres:
  iter+=1
  slope = ( f(x+dx) - f(x-dx) )  / (2*dx)        #find slope
  x = x - damp*(f(x)/slope)
  y=f(x)
  print ('Iteration {} x= {:8.4f}\n'.format(iter,x))

print ('final solution is x= {:8.4f}\n'.format(x))
print ('completed')
