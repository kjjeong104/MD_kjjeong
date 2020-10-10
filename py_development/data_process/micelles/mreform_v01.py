#!/home/kjeong23/softwares/bin/python3.4
# micelle REFORMER. input#1: 1micelle structure input#2: Com point list
# output: resultant gro file
#caution : the 1micelle coordinates should be 'clustered' before operation of this
# algorithm : read 1micelle file info->(clustering is already done)
# -> get COM -> get relative coord from COM -> read COM position file info
# -> locate micelles on each COM -> produce output

import math
import sys
import numpy

def grosplit(str): #own needed function for splitting of .gro format
  splitstr=[str[0:5],str[5:10],str[10:15],str[15:20],str[20:28],str[29:36],str[37:44]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  return splitstr

#Part to load coordinate file, split and put into array
micfile = open(sys.argv[1],'r')
comfile = open(sys.argv[2],'r')
outfile = open(sys.argv[3],'w')
miccrd, comcrd=[],[]
lindex=0
for line in micfile:
  if lindex!=0:
    if lindex==1:
      totnatom=int(line)
      #print(totnatom)
    elif lindex>=2 and lindex<=1+totnatom:
      split=grosplit(line)
      split[0],split[3]=int(split[0]),int(split[3])
      for j in range(4,7):
        split[j]=float(split[j]) #convert strings to their original type
      miccrd.append(split)
    elif lindex==2+totnatom:
      split=line.split()
      micbox=[float(x) for x in split]
  lindex+=1

lindex=0
for line in comfile:
  if lindex!=0:
    if lindex==1:
      totncom=int(line)
    elif lindex>=2 and lindex<=1+totncom:
      split=grosplit(line)
      for j in range(4,7):
        split[j]=float(split[j])
      comcrd.append(split[4:7])
    elif lindex==2+totncom:
      split=line.split()
      combox=[float(x) for x in split]
  lindex+=1

#dictionary for atomic/united atomic masses
amass={'C1':15.035, 'C2':14.027, 'C3': 14.027, 'C4': 14.027, 'C5': 14.027, \
'C6': 14.027, 'C7': 14.027, 'C8': 14.027, 'C9': 14.027, 'CA': 14.027, \
'P1': 30.9738, 'OD1': 15.9994, 'OD2': 15.9994, 'OD3': 15.9994}

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

def calcom(crdentry): 
  #fxn to calculate COM of crd file. Cannot consider pbc. should be pre-clustered
  #masses should be also considered by using dictionary!
  c=[0.0, 0.0, 0.0]
  totmass=0.0
##part 1 : calculating COM
  for row in crdentry:
    m1=amass[row[2]] #corresponding atom's mass
    c=[c[0]+row[4]*m1, c[1]+row[5]*m1, c[2]+row[6]*m1]
    totmass+=m1
  c=[c[0]/totmass,c[1]/totmass,c[2]/totmass]
  return c

#main fxn:main administrator of algorithm progression
# algorithm : read 1micelle file info->(clustering is already done)
# -> get COM -> get relative coord from COM -> read COM position file info
# -> locate micelles on each COM -> produce output
def main():
  com=calcom(miccrd) #get COM of micelle coord
  #transpose miccrd wrt its com 
  for i in range(totnatom):
    for j in range(4,7):
      miccrd[i][j]-=com[j-4]

  #output writing section
  outindex,micindex,nmol=0,0,miccrd[totnatom-1][0]
  outfile.write('Repacked micelles: from coord {}\n'.format(sys.argv[1]))
  outfile.write(' {}\n'.format(totnatom*totncom))
  for comrow in comcrd:
    micindex+=1
    for i in range(totnatom):
      outindex+=1
      x,y,z=comrow[0]+miccrd[i][4],comrow[1]+miccrd[i][5],comrow[2]+miccrd[i][6] 
      outfile.write('{:5}{:5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}\n'.format(miccrd[i][0]+nmol*(micindex-1),miccrd[i][1],miccrd[i][2],outindex,x,y,z))
  outfile.write('   {}   {}   {}'.format(combox[0],combox[1],combox[2]))

if __name__ == "__main__": main()

