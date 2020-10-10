#!/home/kjeong23/softwares/bin/python3.4
# micelle COM closest lattice point assignment program
#<algorithm>
#input file : perfect sigma phase coord file, micelle COM trajectory file(gro of multiple snapshots), micelle morphology text file (**format is from micelle_v051.py)
#(1) From perfect sigma phase coord, make a numpy array of (lattice point index, x, y, z) assignment.
#(2) start reading micelle COM trajectory. For a step : detect number of COMs. If doesn't match with perfect sigma phase, skip that step.
#if matches with perfect sigma phase, read coords/boxsize and store for 1 step. For each COM, assign the index of 'closest perfect sigma phase lattice point'
#(3) read corresponding number of lines (or skip corresponding number of lines) in morphology text file. Rearrange morphology lines, add 'lattice point# ~~' in front of each line.
#now grep command can search for a particular lattice point. In the output file, print new format of morphology description that 'latticepoint index' is added.
#(4) go back to (2) for a new step.

import math
import sys
import numpy

def gro_xyzsplit(str): #own needed function for splitting of .gro format. V05.
  splitstr=[str[20:28],str[28:36],str[36:44]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  for i in range(0,3):
    splitstr[i]=float(splitstr[i])
  return splitstr

def latticeassign(r1,box,lat): #distance calc, find closest one among lattice points
  #lat : 2D array of xyz values of all lattice points
  mindist=box[0]
  minindex=0
  for i in range(lat.shape[0]): #i: index for lattice point
    r=lat[i]-r1
    for j in range(3):
      if r[j]>(box[j]/2.0):
        r[j]-=box[j]
      elif r[j]<(-box[j]/2.0):
        r[j]+=box[j]
    dist=numpy.sqrt(numpy.dot(r,r))
    if dist<mindist:
      mindist=dist
      minindex=i
  minindex+=1 #to make the index starts from 1(not from zero)
  return minindex

#main fxn
def main():
 
  #Load input files
  sigfile = open(sys.argv[1],'r') 
  trjfile = open(sys.argv[2],'r')
  morfile = open(sys.argv[3],'r')
  outfile = open(sys.argv[4],'w') #output file for lattice point assigment

  #construct ideal sigma phase lattice points info array
  sigcrd=numpy.empty((0,3),float)
  lindex=0
  for line in sigfile:
    if lindex !=0:
      if lindex==1:
        totlatp=int(line)
      elif lindex>=2 and lindex<=1+totlatp:
        split=gro_xyzsplit(line)
        sigcrd=numpy.vstack((sigcrd,split))
    lindex+=1

  #start the loop of 'COM trajectory reading'
  #if nmic!=totlatp: skip the whole step
  sindex,lindex=0,0
  actcrd=numpy.empty((0,3),float)
  latlist=[] #lattice point index list
  for line in trjfile:
    if lindex!=0:
      if lindex==1:
        nmic=int(line)
      elif lindex>=2 and lindex<=1+nmic and nmic==totlatp:
        split=gro_xyzsplit(line)
        actcrd=numpy.vstack((actcrd,split))
      elif lindex==2+nmic and nmic==totlatp:
        split=line.split()
        box=numpy.array([float(x) for x in split])

      elif lindex>=3+nmic: #conclusion for 1 step, and initialization again
        if nmic==totlatp:
          morline=morfile.readline()
          outfile.write(morline)
          for row in actcrd:
            latnum=latticeassign(row,box,sigcrd)
            morline=morfile.readline()
            outstr='lattp# ' + str(latnum) + ' ' + morline
            outfile.write(outstr)
          morline=morfile.readline()
        elif nmic!=totlatp:
          for i in range(nmic+2):
            morline=morfile.readline()
          
        #initialzation
        print('step {} complete'.format(sindex))
        actcrd=numpy.empty((0,3),float)
        latlist=[]
        sindex+=1
        lindex=0
    lindex+=1

  sigfile.close()
  trjfile.close()
  morfile.close()
  outfile.close()

if __name__ == "__main__": main()

