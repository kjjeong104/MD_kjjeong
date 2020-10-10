#!/home/kjeong23/softwares/bin/python3.4
# prototype program for micelle detection from coordinate file
# algorithm: get coord&info -> starts loop. -> distance test(consider pbc)
# -> include -> further test again&again -> full 1micelle list -> store
# -> start looking for another micelle
# progress:

import math
import sys
#import numpy

def grosplit(str): #own needed function for splitting of .gro format
  splitstr=[str[0:10],str[10:15],str[15:20],str[20:28],str[29:36],str[37:44]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  return splitstr

#Part to load coordinate file, split and put into array
crdfile = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w') #output file for reduced expression of micelle
minagg=int(input("minimum aggregation number to be recognized as a micelle?recomm:6\n"))
thr=float(input("threshold for micelle detection by tail position? recomm:0.9\n"))
surfstr='DEP'
crd=[]
totcrd=[]
nsurf=0
lindex=0
for line in crdfile:
  if lindex!=0:
    if lindex==1:
      totnatom=int(line)
      print(totnatom)
    elif lindex>=2 and lindex<=1+totnatom:
      split=grosplit(line)
      split[2]=int(split[2])
      for j in range(3,6):
        split[j]=float(split[j]) #convert strings to their original type
      totcrd.append(split)
      if split[1]=='C1':
        nsurf+=1
        crd.append(split)
    elif lindex==2+totnatom:
      split=line.split()
      box=[float(x) for x in split]
  lindex+=1

#dictionary for atomic/united atomic masses
amass={'C1':15.035, 'C2':14.027, 'C3': 14.027, 'C4': 14.027, 'C5': 14.027, \
'C6': 14.027, 'C7': 14.027, 'C8': 14.027, 'C9': 14.027, 'CA': 14.027, \
'P1': 30.9738, 'OD1': 15.9994, 'OD2': 15.9994, 'OD3': 15.9994}

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
  #differently from v02, totcrd of all atoms are considered!
  #masses should be also considered by using dictionary!
  c=[0.0, 0.0, 0.0]
  totmass=0.0
  #refinement of coordinate for wall-crossing micelles
  tends=[] #array for tail end atom indexes
  moveends=[] #array for 'marked molecules' to be moved (negative box direction)
  for i in list:
    if totcrd[i][1]=='C1':
      tends.append(i) #record tail ends. wall-crossing examination will be done
  for i in tends: #wall-crossing test: done for only tail ends.
    for j in tends:
      for k in range(3,6):
        if (totcrd[j][k]-totcrd[i][k])>(box[k-3]/2.0): #if j is cut by (1b)edge
          if [totcrd[j][0],k-3] not in moveends:
            moveends.append([totcrd[j][0],k-3])
        elif (totcrd[i][k]-totcrd[j][k])>(box[k-3]/2.0): #if i is cut by (1b)edge
          if [totcrd[i][0],k-3] not in moveends:
            moveends.append([totcrd[i][0],k-3])

  for x in moveends:
    for y in list:
      if totcrd[y][0]==x[0]:
        totcrd[y][x[1]+3]-=box[x[1]] #negative shift for wall-crossing micelle atoms

  for x in list:
    m1=amass[totcrd[x][1]] #corresponding atom's mass
    c=[c[0]+totcrd[x][3]*m1, c[1]+totcrd[x][4]*m1, c[2]+totcrd[x][5]*m1]
    #print(c[2]+totcrd[x][5]*m1) debug
    totmass+=m1
  c=[round(c[0]/totmass,3),round(c[1]/totmass,3),round(c[2]/totmass,3)]
  print ('calculated center of mass : {} {} {}'.format(c[0],c[1],c[2]))
  return c

#main fxn:main administrator of algorithm progression
def main():
  totmic,already=[],[]
  for i in range(nsurf):#If not contained: starting a new micelle
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
          for j in range(nsurf):
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
  for row in refmic: #tailend member list of 1 micelle
    onemic=[]
    for i in row:
      for j in totcrd:
        if ((str(i+1)+surfstr)==j[0]): #search for specific molID
          onemic.append(j[2]-1)
    c=center(onemic,box)
    coms.append(c)

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
  outfile.write('  {}  {}  {}'.format(box[0],box[1],box[2]))

if __name__ == "__main__": main()

#debug section
#print (natom)
#print (crd)
#print (box)
