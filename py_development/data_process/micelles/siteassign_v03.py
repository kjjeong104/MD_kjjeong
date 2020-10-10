#!/home/kjeong23/softwares/bin/python3.4
# micelle COM closest lattice point assignment program (for sigma phase only!)
#<algorithm>
#input file : perfect sigma phase coord file, micelle COM trajectory file(gro of multiple snapshots), micelle morphology text file (**format is from micelle_v051.py)
#(1) From perfect sigma phase coord, make a numpy array of (lattice point index, x, y, z) assignment.
#(2) start reading micelle COM trajectory. For a step : detect number of COMs. If doesn't match with perfect sigma phase, skip that step.
#if matches with perfect sigma phase, read coords/boxsize and store for 1 step. For each COM, assign the index of 'closest perfect sigma phase lattice point'
#(3) read corresponding number of lines (or skip corresponding number of lines) in morphology text file. Rearrange morphology lines, add 'lattice point# ~~' in front of each line.
#now grep command can search for a particular lattice point. In the output file, print new format of morphology description that 'latticepoint index' is added.
#(4) go back to (2) for a new step.
#v03 : coms file rearrangement : index accordingly (lattice site indices)
#v02 : formatting improvement: lattp# constant indentation, display elapsed time instead of mic#,
# reduce spacing for r= and A= . erase brackets. change right bracket(]) into comma.
#to do that : from morfile, cut position 0~8, 19~20, 31~32  

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
  outfile = open(sys.argv[4],'w') #output file for lattice point assignment
  reafile = open(sys.argv[5],'w') #reindexed com trajectory file

  initt=float(input("Initial time of this trajectory (in ns)?\n"))
  ssize=float(input("Timestep between snapshots (in ns)? ex) 0.04\n"))

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
  reacrd=numpy.zeros((30,3),float) #reindexed array
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
          tns=initt+ssize*sindex #time in ns
          morline=morfile.readline()
          outfile.write(morline)
          for row in actcrd:
            latnum=latticeassign(row,box,sigcrd)
            reacrd[latnum-1]=row
            preline=morfile.readline()
            morline=preline[9:19]+preline[21:31]+preline[33:]
            morline=morline.replace("]",",")
            morline=morline.replace("["," ")
            outstr='lattp# '+"{:2} ".format(latnum)+"{:7.3f} ns ".format(tns)+morline
            outfile.write(outstr)
          morline=morfile.readline()
          #rearranged com file writing
          reafile.write('COmics reindexed: step {} t= {:7.3f} ns\n'.format(sindex,tns))
          reafile.write(' {}\n'.format(nmic))
          outindex=0
          for rearow in reacrd:
            outindex+=1
            reafile.write('{:5}{:5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}\n'.format(outindex,'COM','X',outindex,rearow[0],rearow[1],rearow[2]))
          reafile.write('      {}  {}  {}\n'.format(box[0],box[1],box[2]))
        elif nmic!=totlatp:
          for i in range(nmic+2):
            morline=morfile.readline()
          
        #initialzation
        if (sindex%50)==0:
          print('step {} complete'.format(sindex))
        actcrd=numpy.empty((0,3),float)
        reacrd=numpy.zeros((30,3),float)
        latlist=[]
        sindex+=1
        lindex=0
    lindex+=1

  sigfile.close()
  trjfile.close()
  morfile.close()
  outfile.close()
  reafile.close()

if __name__ == "__main__": main()

