#!/home/kjeong23/softwares/bin/python3.4

import math
import numpy
import numpy.linalg
from numpy.linalg import inv
# Pulay mixing method for numerical solution search practice
# Basically rely on Least-square minimization of error

#Variables section
thres=0.0001  #threshold
target=0      #target value for the function
guess=0     #initial guess clue for x
dx=0.1    #increment of x to measure slope of f(x). search for both direction
nhist=5   # number of history for pulay mixing
#damp=0.7  #damping factor

#solve equation f(x)=c, when c is the target.
# for this case, let f(x)=x^3-2x+2

def f(x):
  value=x**3-2*x+2
  return value

def matrixB(vec):
  B = numpy.zeros((vec.shape[0]+1,vec.shape[0]+1))
  for i in range(vec.shape[0]+1):
    for j in range(vec.shape[0]+1):
      if i<vec.shape[0] and j<vec.shape[0]:
        B[i][j]=vec[i]*vec[j]
      elif i==vec.shape[0] and j==vec.shape[0]:
        B[i][j]=0
      else:
        B[i][j]=-1

  print(B)
  return B

x=guess
y=f(x)
iter=0

while abs(y-target) > thres:
  iter+=1
  if iter<=1:
    xvec=numpy.array([x])    #set the vector of historical x. this will be appended
    slope = ( f(x+dx) - f(x-dx) )  / (2*dx)        #find slope
    x = x - (f(x)/slope)
    xvec=numpy.append(xvec,x) #now xvec contains x1, x2. Can't do mixing until here
  else:
    if iter<nhist:
      slope = ( f(x+dx) - f(x-dx) )  / (2*dx)
      newx= x - (f(x)/slope)  #candidate of new x before applying mixing
      dxvec=numpy.zeros(iter) #delta x array. this is included in least square problem
      for i in range(xvec.shape[0]):
        if i<xvec.shape[0]-1:
          dxvec[i]=xvec[i+1]-xvec[i]
        else:
          dxvec[i]=newx-xvec[i]

      #construction of dxvec is complete
      B = matrixB(dxvec)
      cvec=numpy.zeros(B.shape[0])
      cvec[B.shape[0]-1]=-1

      Binv= inv(B)
      avec=Binv.dot(cvec) # coefficient vector

      avec=numpy.delete(avec,avec.shape[0]-1)
      x=numpy.dot(avec,xvec) #the "true" new x
      xvec=numpy.append(xvec,x)

    else:
      slope = ( f(x+dx) - f(x-dx) )  / (2*dx)
      newx= x - (f(x)/slope)  #candidate of new x before applying mixing
      dxvec=numpy.zeros(nhist) #delta x array. this is included in least square problem
      for i in range(nhist):
        if i<nhist-1:
          dxvec[i]=xvec[i+1]-xvec[i]
        else:
          dxvec[i]=newx-xvec[i]

      B = matrixB(dxvec)
      cvec=numpy.zeros(B.shape[0])
      cvec[B.shape[0]-1]=-1

      Binv= inv(B)
      print(Binv)
      avec=Binv.dot(cvec) # coefficient vector

      avec=numpy.delete(avec,avec.shape[0]-1)
      print('avec is {} \n'.format(avec))
      x=numpy.dot(avec,xvec) #the "true" new x
      #push back xvec elements
      xvec=numpy.append(xvec,x)
      xvec=numpy.delete(xvec,0)

  y=f(x)
  print ('Iteration {} x= {:8.4f}\n'.format(iter,x))

print ('final solution is x= {:8.4f}\n'.format(x))
print ('completed')
