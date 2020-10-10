#!/home/kjeong23/softwares/bin/python3.4
# micelle COM tracker program (good for self-assembled system)
#<algorithm>
#input file : micelle COM trajectory file(gro of multiple snapshots), micelle morphology text file (**format is from micelle_v051.py)
#(1) First user puts designated number of micelles want to set for tracking.
#(2) Start reading micelle COM trajectory. When it detects the step that # micelles meet the requirement,
# record the coordinates,box. If #COMs doeesn't match with requirement, skip that step.
#(3) For each COM, assign the index of 'most similar member micelle'. (important change of address assigning algorithm!)
#(4) read corresponding number of lines (or skip corresponding number of lines) in morphology text file. Rearrange morphology lines, add 'lattice point# ~~' in front of each line. 

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

'''
def latticeassign(r1,box,lat): #distance calc, find closest one among lattice points
  lat : 2D array of xyz values of all lattice points
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
'''

def memtrack(target,fullentry) #by comparison of members, find correct address of a micelle
  #criterion : highest similarity
  #target: memberlist of 1 micelle, fullentry : memberlist of all micelles (from previous step)
  maxagree,index=0,0
  for i in range(len(fullentry)):
    agree=0
    for x in target:
      if x in fullentry[i]:
        agree+=1
    if agree>maxagree
      maxagree=agree
      index=i
  index+=1 #to make index starts from 1 (not from 0)
  return index

#main fxn
def main():
 
  #Load input files
  comfile = open(sys.argv[1],'r')
  morfile = open(sys.argv[2],'r')
  outfile = open(sys.argv[3],'w') #output file for lattice point assignment
  reafile = open(sys.argv[4],'w') #reindexed com trajectory file

  initt=float(input("Initial time of this trajectory (in ns)?\n"))
  ssize=float(input("Timestep between snapshots (in ns)? ex) 0.04\n"))
  totlatp=int(input("Designated # of micelles for tracking? ex) 33 \n"))

  #start the loop of 'COM trajectory reading'
  #if nmic!=totlatp: skip the whole step
  sindex,lindex,flag=0,0,0
  actcrd=numpy.empty((0,3),float)
  reacrd=numpy.zeros((totlatp,3),float) #reindexed array
  fullentry=[] #lattice point index list
  for line in comfile:
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
          if flag==0: #first suitable step
            tns=initt+ssize*sindex #time in ns
            morline=morfile.readline()
            outfile.write(morline)
            flag,i=1,0
            for row in actcrd:
              preline=morfile.readline()
              morline=preline[9:19]+preline[21:31]+preline[33:]
              morline=morline.replace("]",",")
              morline=morline.replace("["," ")
              memline=preline[84:]
              memline=memline.replace("]","")
              memline=memline.replace("[","")
              memline=memline.replace(",","")
              onemmem=memline.split()
              fullentry.append(onemmem)
              reacrd[i]=row
              outstr='lattp# '+"{:2} ".format(i+1)+"{:7.3f} ns ".format(tns)+morline
              outfile.write(outstr)
              i+=1
            morline=morfile.readline()
          elif: #not the first step
            tns=initt+ssize*sindex #time in ns
            morline=morfile.readline()
            outfile.write(morline)
            newentry=fullentry
            for row in actcrd:
              preline=morfile.readline()
              morline=preline[9:19]+preline[21:31]+preline[33:]
              morline=morline.replace("]",",")
              morline=morline.replace("["," ")
              memline=preline[84:]
              memline=memline.replace("]","")
              memline=memline.replace("[","")
              memline=memline.replace(",","")
              target=memline.split()
              latnum=memtrack(target,fullentry)
              newentry[latnum-1]=target
              reacrd[latnum-1]=row
              outstr='lattp# '+"{:2} ".format(latnum)+"{:7.3f} ns ".format(tns)+morline
              outfile.write(outstr)
            morline=morfile.readline()
            fullentry=newentry #renew the fullentry
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
        if (sindex%20)==0:
          print('step {} complete'.format(sindex))
        actcrd=numpy.empty((0,3),float)
        reacrd=numpy.zeros((totlatp,3),float)
        sindex+=1
        lindex=0
    lindex+=1

  comfile.close()
  morfile.close()
  outfile.close()
  reafile.close()

if __name__ == "__main__": main()

