#!/home/kjeong23/softwares/bin/python3.4
# prototype program for micelle detection from coordinate file
# algorithm: get coord&info -> starts loop. -> distance test(consider pbc)
# -> include -> further test again&again -> full 1micelle list -> store
# -> start looking for another micelle
# progress:

import math
import sys
#import numpy

#Part to load coordinate file, split and put into array
crdfile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w') #output file for reduced expression of micelle
minagg=int(input("minimum aggregation number to be recognized as a micelle?recomm:6\n"))
thr=float(input("threshold for micelle detection by tail position? recomm:0.9\n"))
crd=[]
for line in crdfile:
  crd.append(line.split())
natom=crd.pop([1][0])
natom=int(natom[0])
del crd[0]
box=crd.pop()
box=[float(x) for x in box]
for i in range(natom):
  crd[i][2]=int(crd[i][2])
  for j in range(3,6):
    crd[i][j]=float(crd[i][j])  #convert strings to their original type

def distdet(r1,r2,thr,box): #distance checking function. returns true for detection
  r=[r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]]
  for i in range(3):
    if r[i]>(box[i]/2.0):
      r[i]-=box[i]
    elif r[i]<(-box[i]/2.0):
      r[i]+=box[i]
  dist=math.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
  if dist<=thr:
    det=True 
  else:
    det=False
  return det

def lookup1d(n,array):#1d array lookup for a particular number
  det=False
  if n in array:
    det=True
  return det

def lookup2d(n,array):#2d array lookup for a particular number
  det=False
  for row in array:
    if n in row:
      det=True
  return det

def center(list,box): #fxn to calculate COM of an atom group, considering pbc.
  # only crd of tailends are considered
  c=[0.0, 0.0, 0.0]
  #refinement of coordinate for wall-crossing micelles
  for i in list:
    for j in list:
      for k in range(3,6):
        if (crd[j][k]-crd[i][k])>(box[k-3]/2.0):
          crd[j][k]-=box[k-3]
        elif (crd[i][k]-crd[j][k])>(box[k-3]/2.0):
          crd[i][k]-=box[k-3]

  for x in list:
    c=[c[0]+crd[x][3], c[1]+crd[x][4], c[2]+crd[x][5]]
  n=len(list)
  c=[round(c[0]/n,3),round(c[1]/n,3),round(c[2]/n,3)]
  print ('calculated center of mass : {} {} {}'.format(c[0],c[1],c[2]))
  return c

#main fxn:main administrator of algorithm progression
def main():
  totmic,already=[],[]
  for i in range(natom):#If not contained: starting a new micelle
    onemic=[]
    if lookup2d(i,totmic)==False:
      onemic.append(i)
      search=1
      print('creating new micelle,molindex {}'.format(i)) ####debug
    while search==1:#starting of repetitive searching algorithm
      search=0
      for k in onemic:#stops when increase of number in onemic stopped
        if lookup1d(k,already)==False:
          #print('having primary target,molindex {}'.format(k)) ####debug
          for j in range(natom):
            if j!=k and lookup2d(j,totmic)==False and lookup1d(j,onemic)==False:
              if(distdet(crd[k][3:6],crd[j][3:6],thr,box)): #detection
                onemic.append(j)
                search=1
          already.append(k)
    if len(onemic)!=0:
      totmic.append(onemic)
    already=[]

  #list refinement section
  refmic=[]
  hom=[]
  unagg=0
  for row in totmic:
    if len(row)>=minagg:
      refmic.append(row)
    else:
      hom.append(row)
      unagg+=len(row)

  #com reduction
  coms=[]
  for row in refmic:
    c=center(row,box)
    coms.append(c)
  #coms merging for erroneous detection
  black=[] #blacklist micelles
  k=len(coms)
  for i in range(k):
    for j in range(i+1,k):
      #check if two coms are too close to each other
      if(distdet(coms[i],coms[j],thr*1.5,box)): 
        #if lookup1d(i,black)==False:
        black.append(i)
        #if lookup1d(j,black)==Flase:
        black.append(j)
        c=center(refmic[i]+refmic[j],box)
        coms.append(c)
        refmic.append(refmic[i]+refmic[j])
  #merging:last step
  black.sort(reverse=True)
  for i in black:
    refmic[i]='!'
    coms[i]='!'
  for i in range(0,refmic.count('!')):
    refmic.remove('!')
  for i in range(0,coms.count('!')):
    coms.remove('!')

  #printing section
  micindex=0
  for row in refmic:
    micindex+=1
    print ('micelle# {} aggnum {} members: {}'.format(micindex,len(row),row))
  print ('{} Unaggregated molecules: {}'.format(unagg,hom))

  #com output writing section
  outindex=0
  outfile.write('Centers of micelles: from coord {}\n'.format(sys.argv[1]))
  outfile.write(' {}\n'.format(micindex))
  for row in coms:
    outindex+=1
    outfile.write('{:5}{:5}{:5}{:5}{:8}{:8}{:8}\n'.format(outindex,'COM','X',outindex,row[0],row[1],row[2]))
  outfile.write(' {} {} {}'.format(box[0],box[1],box[2]))

if __name__ == "__main__": main()

#debug section
#print (natom)
#print (crd)
#print (box)
