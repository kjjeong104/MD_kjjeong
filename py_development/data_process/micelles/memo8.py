#!/home/kjeong23/softwares/bin/python3.4

import math
import numpy
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
#fcrd=numpy.array([[0,1,0],[4,0,0],[0,4,0]])
#apex=numpy.array([1,8,6])

#r12,r23,r31=fcrd[0]-fcrd[1],fcrd[1]-fcrd[2],fcrd[2]-fcrd[0]
#sa,sb,sc=numpy.linalg.norm(r12),numpy.linalg.norm(r23),numpy.linalg.norm(r31)
#ss=(sa+sb+sc)/2.0
#oneS=numpy.sqrt( ss * (ss-sa) * (ss-sb) * (ss-sc) )

#Volume of trigonal pyramid is (1/3)*A*H.
#pnorm=numpy.cross(r23,r12)
#pnorm=pnorm/( numpy.linalg.norm(pnorm) )
#r14=fcrd[0]-apex
#h=numpy.linalg.norm( numpy.dot(pnorm,r14) ) /numpy.linalg.norm(pnorm)
#oneV=oneS*h/3.0

#SandV=numpy.array([oneS,oneV])
#print (pnorm)
#print (SandV)

miccrd=numpy.array([ [0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,0,1],[1,0,1],[0,1,1],[1,1,1] ])
#tri=Delaunay(cubecrd)
#hull=ConvexHull(cubecrd)

#print(tri.simplices)
#print(hull.simplices)

def Tpyramid(fcrd,apex): #triangle area and trigonal pyramid volume calculator.
  #fcrd:3*3 array of triangle points coords. apex:1*3 array of apex coord. should be pbc refined.
  #Get area of triangle by Heron's formula T=sqrt(s*(s-a)*(s-b)*(s-c)) when s=(a+b+c)/2
  #a=r12, b=r23, c=r31.
  r12,r23,r31=fcrd[0]-fcrd[1],fcrd[1]-fcrd[2],fcrd[2]-fcrd[0]
  sa,sb,sc=numpy.linalg.norm(r12),numpy.linalg.norm(r23),numpy.linalg.norm(r31)
  ss=(sa+sb+sc)/2.0
  oneS=numpy.sqrt( ss * (ss-sa) * (ss-sb) * (ss-sc) )
  #Volume of trigonal pyramid is (1/3)*A*H.
  pnorm=numpy.cross(r12,r23)
  #pnorm=pnorm/( numpy.linalg.norm(pnorm) )
  r14=fcrd[0]-apex
  h=numpy.linalg.norm( numpy.dot(pnorm,r14) ) /numpy.linalg.norm(pnorm)
  oneV=oneS*h/3.0
  SandV=numpy.array([oneS,oneV])
  return SandV

def proxtri(crd): #proximity triangle set constructor
  tcrd=numpy.empty((0,3))
  hull=ConvexHull(crd)
  hullsim=hull.simplices

  plt.triplot(crd[:,0], crd[:,1], hull.simplices.copy())
  plt.plot(crd[:,0], crd[:,1], 'o')
  plt.show()

  for row in hullsim:
    for x in row:
      tcrd=numpy.vstack((tcrd,crd[x]))

  return tcrd

#Get geometric midpoint (this will serve as common apex of triangular pyramids)
mid=numpy.mean(miccrd,axis=0)

#Construct proximity triangles, and make them as 3N_triangle*3 array.
#Each 3 rows represents 1 triangle, made of 3 points. (x1,y1,z1),(x2,y2,z2),(x3,y3,z3)
tricrd=proxtri(miccrd)

S,V=0,0
for t in range(int(tricrd.shape[0]/3)): #t: "triangle index"
  currt,nextt=3*t,3*t+3 #row index for current triangle, next triangle
  oneSV=Tpyramid(tricrd[currt:nextt],mid)
  S+=oneSV[0]
  V+=oneSV[1]

Q=36.0 * math.pi *V*V/(S*S*S)
micinfo=numpy.array([S,V,Q])
print(micinfo)

