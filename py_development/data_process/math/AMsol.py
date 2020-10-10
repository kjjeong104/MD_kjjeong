#!/home/kjeong23/softwares/bin/python3.4

import math
import numpy
import numpy.linalg
from numpy.linalg import inv
from scipy.linalg import solve
# Anderson mixing method for numerical solution search practice(DIIS)
# Basically rely on Least-square minimization of error

#Variables section
thres=0.0001  #threshold
target=0      #target value for the function
guess=0     #initial guess clue for x
dx=0.05    #increment of x to measure slope of f(x). search for both direction
nhist=5   # number of history for pulay mixing
#damp=0.7  #damping factor

#solve equation f(x)=c, when c is the target.
# for this case, let f(x)=x^3-2x+2

def f(x):
  value=x**3-2*x+2
  return value

def matrixU(dvec,s):
  #s : last index of the error vector dvec, size of U
  U = numpy.zeros((s,s))
  for i in range(s):
    for j in range(i,s):
      U[i][j]=(dvec[s]-dvec[s-i-1])*(dvec[s]-dvec[s-j-1])
      U[j][i]=U[i][j]

  print(U)
  return U

def vectorV(dvec,s):
  #s : last index of the error vector dvec, size of V
  V = numpy.zeros(s)
  for i in range(s):
    V[i]=(dvec[s]-dvec[s-i-1])*dvec[s]

  return V

x=guess
iter=0
xhist=numpy.zeros(nhist+1)
dvec=numpy.zeros(nhist+1)

#first, need to run twice to get enough history to apply AM scheme
for i in range(3):
  iter+=1
  if iter>1: #except first step
    x = x-(f(x)/slope) #calc updated x
  slope = ( f(x+dx) - f(x-dx) )  / (2*dx)        #find slope
  d = -(f(x)/slope) #planned NR displacement
  y = f(x)
  #d = y - target #deviation
  xhist[iter-1]=x
  dvec[iter-1]=d
  print ('Iteration {} x= {:8.4f}\n'.format(iter,x))

iter-=1 #correct iter index for the other loop

while abs(y-target) > thres:
  iter+=1
  xx,dd=0,0

  if iter<=nhist:
    s = iter-1
  else:
    s = nhist

  #Do the matrix algebra part
  U = matrixU(dvec,s)
  V = vectorV(dvec,s)
  #use SVD to solve eqsystem Uc=v.
  #UM,DM,VM = numpy.linalg.svd(U)
  #Uinv= numpy.dot(VM.transpose(), numpy.dot(numpy.diag(DM**-1),UM.transpose()))
  #cvec=Uinv.dot(V)
  cvec=solve(U,V)
  print(cvec)

  xx=x
  dd=dvec[s]
  for i in range(s):
    xx+=cvec[i]*(xhist[s-i-1]-xhist[s])
    dd+=cvec[i]*(dvec[s-i-1]-dvec[s])

  x=xx+dd
  y=f(x)

  #dev filling
  if iter<=nhist:
    slope = ( f(x+dx) - f(x-dx) )  / (2*dx)
    print(f(x+dx), f(x-dx), 2*dx)
    dvec[iter]= -(f(x)/slope)
    xhist[iter]=x
  else:
    #pushing
    xhist[0:nhist-1]=xhist[1:nhist]
    dvec[0:nhist-1]=dvec[1:nhist]

    slope = ( f(x+dx) - f(x-dx) )  / (2*dx)
    dvec[nhist]= -(f(x)/slope)
    xhist[nhist]=x

  print ('Iteration {} x= {:8.4f}\n'.format(iter,x))

print ('final solution is x= {:8.4f}\n'.format(x))
print ('completed')

