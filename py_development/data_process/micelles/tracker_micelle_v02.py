#!/home/kjeong23/softwares/bin/python3.4
# micelle COM tracker program (good for self-assembled system)
#<algorithm>
#input file : micelle COM trajectory file(gro of multiple snapshots), micelle morphology text file (**format is from micelle_v051.py)
#(1) First user puts designated number of micelles want to set for tracking.
#(2) Start reading micelle COM trajectory. When it detects the step that # micelles meet the requirement,
# record the coordinates,box. If #COMs doeesn't match with requirement, skip that step.
#(3) For each COM, assign the index of 'most similar member micelle'. (important change of address assigning algorithm!)
#(4) read corresponding number of lines (or skip corresponding number of lines) in morphology text file. Rearrange morphology lines, add 'lattice point# ~~' in front of each line.
#v02 : simplify input (only COM coord file). consider dislocation of micelles. for micelle fission: preserve unaffected micelles.
#code abandoned. Decided to embed micelle tracking code in kinmic_v04 series.

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

def pbcdist(r1,r2,box)
  dim=box.shape[0]
  r=r2-r1
  for k in range(dim):
    if r[k]>(box[k]/2.0):
      r[k]-=box[k]
    elif r[k]<(box[k]/2.0):
      r[k]+=box[k]
  dist=numpy.sqrt(numpy.dot(r,r))
  return dist

def tracker(prevcrd,stepcrd,box,oldiline)
  #first, simplest tracking
  nmic1,nmic2=prevcrd.shape[0],stepcrd.shape[0]
  flag=1
  movlist=[]
  if nmic1==nmic2:
    for i in range(nmic1):
      if pbcdist(prevcrd[i],stepcrd[i],box)>=0.85: movlist.append(i)
    if len(movlist)==0: #no significant translocation
      matchline=list(oldiline)
    else: #if swapping happened, rigorous search loop
      mindist=box[0]
      swapresult=[]
      for i in movlist:
        for j in range(nmic1):
          dist=pbcdist(prevcrd[i],stepcrd[j],box)
          if dist<=mindist:
            mindist=dist
            minindex=j
        swapresult.append(minindex)
      matchline=list(oldiline)
      for i in movlist:
        matchline[movlist[i]]=oldiline[swapresult[i]]
  else: #if fission/fusion exists, rigorous search loop
    mindist=box[0]
    swapresult=[]
    for i in range(nmic1):
      for j in range(nmic2):
        dist=pbcdist(prevcrd[i],stepcrd[j],box)
        if dist<=mindist:
          mindist=dist
          minindex=j
      swapresult.append(minindex)
  return matchline

#main fxn
def main():
 
  #Load input files
  comfile = open(sys.argv[1],'r') #micelle v061 above version COM output
  outfile = open(sys.argv[2],'w') 
  #output v02: tracking movement of micelles, print "micelle index list" in the order of matching result with previous snapshots

  #start the loop of 'COM trajectory reading'
  sindex,lindex,flag=0,0,0
  stepcrd=numpy.empty((0,3),float)
  for line in comfile:
    lindex+=1
    if lindex==2:
      nmic=int(line)
    elif lindex>=3 and lindex<=2+nmic:
      split=gro_xyzsplit(line)
      stepcrd=numpy.vstack((stepcrd,split))
    elif lindex==3+nmic: #last line of 1 step
      split=line.split()
      box=numpy.array([float(x) for x in split])
      if sindex==0: #first step
        prevcrd=numpy.copy(stepcrd)
        oldiline=[int(x+1) for x in range(nmic)]
      matchline=tracker(prevcrd,stepcrd,box,oldiline)
      sindex+=1
      lindex=0

  comfile.close()
  outfile.close()

if __name__ == "__main__": main()

