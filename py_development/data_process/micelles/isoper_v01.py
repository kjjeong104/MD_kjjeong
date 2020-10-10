#!/home/kjeong23/softwares/bin/python3.4
# Isoperimetric quotient calculation program (for phosphonate micellar system)
# algorithm: get atomindex info from grofile -> read ASCII traj(1step) and morfile(1step) together
# -> set list of P atoms for each micelle -> regarding pos, construct 'triangles' to enclose micelle.
# -> add up areas of triangles to get surface area S of micelle.
# -> Then for each triangular face,construct cones with an interior point and sum volume of pyramids.
# -> Calculate Q(isoperimetric quotient) for each lattice site. calc time-average too.
# Isoperimetric quotient Q=36*pi*V^2/(S^3). This is 1 for perfect sphere. Lowers as S increases.
# 4 arguments : grofile ASCIItraj morfile output
# caution : for morfile, use reindexed files to clarify the lattice site of micelle.

import math
import sys
import numpy
from scipy.spatial import ConvexHull
import timeit

def gro_atomtype(str): #own needed function for splitting of .gro format.
  splitstr=str[10:15]
  splitstr=splitstr.replace(" ","")
  return splitstr

def ascii_trajsplit(str): #reading coord from traj ASCII. box size cannot be read by this.
  splitstr=[str[16:28],str[29:42],str[43:56]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  posdata=numpy.array([float(x) for x in splitstr])
  return posdata

def pbc_refine(crd,box): #PBC refining function for wall-crossing clusters
  moveends=[] #marked particles for moving (-box)
  npart=crd.shape[0]
  for i in range(npart):
    for j in range(i,npart):
      for k in range(3):
        if ( crd[j][k] - crd[i][k]  ) > (box[k]/2.0): #if j is cut by (1b)edge
          if [j,k] not in moveends: 
            moveends.append([j,k])
        elif ( crd[i][k] - crd[j][k] ) > (box[k]/2.0): #if i is cut by (1b)edge
          if [i,k] not in moveends:
            moveends.append([i,k])

  for x in moveends:
    crd[x[0]][x[1]]-=box[x[1]] #negative shift for wall-crossing micelle atoms

  return crd

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

  for row in hullsim:
    for x in row:
      tcrd=numpy.vstack((tcrd,crd[x]))

  return tcrd

def isocalc(miccrd,box): #isoperimetric quotient calculator by the method of proximity triangle.
  #first, refine the surface atoms' coordinates considering PBC. They should cluster.
  refcrd=pbc_refine(miccrd,box)
  #Get geometric midpoint (this will serve as common apex of triangular pyramids)
  mid=numpy.mean(refcrd,axis=0)

  #Construct proximity triangles, and make them as 3N_triangle*3 array.
  #Each 3 rows represents 1 triangle, made of 3 points. (x1,y1,z1),(x2,y2,z2),(x3,y3,z3)
  tricrd=proxtri(refcrd)

  S,V=0,0
  for t in range(int(tricrd.shape[0]/3)): #t: "triangle index"
    currt,nextt=3*t,3*t+3 #row index for current triangle, next triangle
    oneSV=Tpyramid(tricrd[currt:nextt],mid)
    S+=oneSV[0]
    V+=oneSV[1]

  Q=36.0 * math.pi *V*V/(S*S*S)
  micinfo=numpy.array([S,V,Q])
  return micinfo

def main():
  #Load files
  grofile = open(sys.argv[1],'r')
  trjfile = open(sys.argv[2],'r')
  morfile = open(sys.argv[3],'r') #reindexed morphology summary file (including members)
  outfile = open(sys.argv[4],'w') #output file for number of isoper- quotient calc for latt sites

  initt=float(input("Initial time of this trajectory (in ns)?\n"))
  fint=float(input("Final time of the whole trajectory in ns? ex) 499.96 \n"))
  ssize=float(input("Timestep between snapshots (in ns)? ex) 0.04\n"))
  nmic=30
 
  start_time=timeit.default_timer()
  #load grofile -> get atomindex list of desired atom first (phosphorus?)
  #be careful: this code only works when only 1 representative atom matches with 1 molecule
  ailist=[]
  lindex,nsurf=0,0
  for line in grofile:
    if lindex !=0:
      if lindex==1:
        totnatom=int(line)
      elif lindex>=2 and lindex<=1+totnatom:
        split=gro_atomtype(line)
        if split=='P1': #track phosphorus atom. ailist index:molecule index(from0),save atomi(from0)
          ailist.append(lindex-2) #ailist consists of numbers, not texts
          nsurf+=1
    lindex+=1

  #prepare bins
  nstep=int((fint-initt)/ssize)
  box=numpy.zeros((nstep+1,3))

  #load ASCII traj, store surface atom coord (phosphorus?)
  #morfile is read later.
  lindex,sindex,mindex=0,0,0
  sacrd=numpy.zeros((nstep+1,nsurf,3)) #surface atom coord
  for line in trjfile:
    if lindex==3: #box vector x
      splitstr=[line[18:30]]
      for i in range(len(splitstr)):
        splitstr[i]=splitstr[i].replace(" ","")
      box[sindex][0]=float(splitstr[0])
    elif lindex==4: #box vector y
      splitstr=[line[31:44]]
      for i in range(len(splitstr)):
        splitstr[i]=splitstr[i].replace(" ","")
      box[sindex][1]=float(splitstr[0])
    elif lindex==5: #box vector z
      splitstr=[line[45:58]]
      for i in range(len(splitstr)):
        splitstr[i]=splitstr[i].replace(" ","")
      box[sindex][2]=float(splitstr[0])
    elif lindex>=7 and lindex<=6+totnatom:
      if (lindex-7) in ailist: #if reading phosphorus coord (lindex-7=aindex)
        split=ascii_trajsplit(line)
        sacrd[sindex][mindex]=split
        mindex+=1
    if lindex==6+totnatom: #conclusion for 1 step, initialize
      mindex=0
      sindex+=1
      lindex=-1
    lindex+=1
  elapsed=timeit.default_timer() - start_time
  print('reading traj step {} time {:8.4f}'.format(sindex,elapsed))

  #morfile processing. Should be careful to check timestep.
  isoinfo,isoave=numpy.zeros((nmic,3)),numpy.zeros((nmic,3)) #info summary(1step,time ave).each row:S,V,Q. row index=mic site#-1
  count,nefstep=0,0 #pristine count of micelles in whole traj. effective # of steps.
  for line in morfile:
    #micelle loop
    if len(line)>=90: #for the line actually telling members
      memline=line[91:]
      memline=memline.replace(",","")
      memsplit=memline.split() #1d list of member surf indices
      memlist=[]
      for x in memsplit: #should change datatype of memsplit into int
        memlist.append(int(x))
      ltsplit=line[6:17].split() #lattice point index and time
      ltsplit[0],ltsplit[1]=int(ltsplit[0]),float(ltsplit[1])
      sindex=int((ltsplit[1]-initt)/ssize) #monitor step index
      #micelle surface atoms' coord stacking
      miccrd=numpy.empty((0,3),float)
      for a in memlist:
        miccrd=numpy.vstack((miccrd,sacrd[sindex][a]))
      #1 micelle surface atoms' crd saved. now calc S,V.
      isoinfo[ltsplit[0]-1]=isocalc(miccrd,box[sindex]) #send micelle surafce crds,boxinfo into Q calculator
      count+=1

      #when a step completes
      if count==nmic:
        count=0
        nefstep+=1
        isoave+=isoinfo
        #1step printing
        outfile.write('step# {} micsite# S V Q\n'.format(sindex))
        for i in range(nmic):
          outfile.write('{:2} {:8.4f} {:8.4f} {:8.4f}\n'.format(i+1,isoinfo[i][0],isoinfo[i][1],isoinfo[i][2]))
        elapsed=timeit.default_timer() - start_time
        print('calculating S,V complete: step# {} time {:8.4f}'.format(sindex,elapsed))

  #final printing section for time average
  #unit of S,V: nm^2,nm^3. Q is unitless.
  isoave=isoave/nefstep
  outfile.write('Final time average over {} effective steps: micsite# S(nm^2) V(nm^3) Q\n'.format(nefstep))
  outfile.write('Caution: S,V,Q Averaging done separately, so S_av, V_av do not reproduce Q_av.\n')
  for i in range(nmic):
    outfile.write('{:2} {:8.4f} {:8.4f} {:8.4f}\n'.format(i+1,isoave[i][0],isoave[i][1],isoave[i][2]))
  
  elapsed=timeit.default_timer() - start_time
  print('averaging S,V complete: time {:8.4f}'.format(elapsed))

  grofile.close()
  trjfile.close()
  morfile.close()
  outfile.close()

if __name__ == "__main__": main()
