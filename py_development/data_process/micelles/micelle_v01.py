#!/home/kjeong23/softwares/bin/python3.4
# prototype program for micelle detection from coordinate file
# algorithm: get coord&info -> starts loop. -> distance test(consider pbc)
# -> include -> further test again&again -> full 1micelle list -> store
# -> start looking for another micelle
# progress:

import math
#import numpy

#Part to load coordinate file, split and put into array
#crdfile = open('TMAC10PO3_rev_w30_200ns_tailend.gro','r')
crdfile = open('TMAy2_C10PO3_40wt_tailend.gro','r')
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
  for j in range(3,9):
    crd[i][j]=float(crd[i][j])  #convert strings to their original type
thr=0.7000 #threshold for micelle detection by tail position
minagg=6 #minimum aggregation number to be recognized as a micelle

def distdet(r1,r2,thr,box): #distance checking function. returns true for detection
  r=[r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]]
  for i in range(3):
    if r[i]>(box[i]/2.0):
      r[i]-=box[i]
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
  for row in totmic:
    if len(row)>=minagg:
      refmic.append(row)
    else:
      hom.append(row)

  #printing section
  micindex=1
  for row in refmic:
    print ('micelle# {} aggnum {} members: {}'.format(micindex,len(row),row))
    micindex+=1
  print ('Unaggregated molecules: {}'.format(hom))

if __name__ == "__main__": main()

#debug section
#print (natom)
#print (crd)
#print (box)
