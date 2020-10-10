#!/home/kjeong23/softwares/bin/python3.4
# micelle kinetics analysis program
# input file : merged lattice point assignment summary file for entire time period
#v02: site-specific kinetic info, micelle exchange portion, exchange rate, gain&lose
#v03: retention time histogram, raw data record for (m-index,psite#,site#,ret.time)
#v03 renovation: also counted diffusion out to solvent media.
#v04:now should be able to count even if clustering data has dislocation/fission.
#also calculate "time dependent relaxation function" of member retention
#policy for fission/fusion : just "uncounting" the broken/coalesced micelles
#** tracker matrix buggy. abandon code ***

import math
import sys
import numpy
import timeit

mthres=6 #threshold aggregation number of an independent micelle : 6
def comparison(a1,a2): #getting retention fraction of micelle during change a1->a2
  coin=0
  for x in a1:
    if x in a2:
      coin+=1
  retf=coin/len(a1)
  return retf

def tracker(prevmor,stepmor,oldiline): #two steps' morphol matrix, prev step's match index list
  #first, simplest tracking
  nmic1,nmic2=len(prevmor),len(stepmor)
  movlist,swapresult=[],[]
  if nmic1==nmic2: #equal # of micelles
    for j in range(nmic2):
      if comparison(prevmor[j],stepmor[j])<0.50 or abs(len(prevmor[j])-len(stepmor[j]))>=mthres: 
        movlist.append(j) #dislocation judgement
    if len(movlist)==0: #no site translocation
      matchline=list(oldiline)
    else: #if translocation happened, rigorous search loop
      for j in movlist:
        maxfract=0
        for i in range(nmic1):
          retf=comparison(prevmor[i],stepmor[j])
          if retf>=maxfract:
            maxfract=retf
            maxindex=i
        swapresult.append(maxindex)
      matchline=list(oldiline)
      for i in range(len(movlist)):
        if abs(len(prevmor[movlist[i]])-len(stepmor[swapresult[i]]))>=mthres: matchline[swapresult[i]]='X'
        else: matchline[swapresult[i]]=oldiline[movlist[i]]
  else: #if fission/fusion happened, rigorous search loop
    #we need to distinguish fissed/fused micelles, and which didn't.
    #include "aggregation number" for the matching criterion
    #loop order is reversed. focus on the latter step.
    #print('prevmor stepmor {} {}, unequal mic tracker {}'.format(prevmor,stepmor,oldiline))
    matchline=[]
    for j in range(nmic2):
      maxfract=0
      for i in range(nmic1):
        retf=comparison(stepmor[j],prevmor[i])
        if retf>=maxfract:
          maxfract=retf
          maxindex=i
      #if the micelle matches the best has far different aggN, that means fission/fusion
      if abs(len(prevmor[maxindex])-len(stepmor[j]))>=mthres: matchline.append('X')
      else: matchline.append(oldiline[maxindex])
  return matchline

def relaxcorr_1step(steplist1,steplist2,tracker1,tracker2): 
  #two snapshots: relaxation fraction function, micelle count
  #finally, relaxation fraction should be "number-averaged" by number of micelles
  #first: try direct mapping of lattice site, if tracker entry is the same.
  nmic1,nmic2=len(steplist1),len(steplist2)
  totrf=0
  if tracker1==tracker2:
    mcount=nmic1
    for i in range(nmic1):
      dmf=comparison(steplist1[i],steplist2[i]) #directly matched fraction
      totrf+=dmf
  else:
    mcount=0
    for i in range(nmic1):
      for j in range(nmic2):
        if tracker1[i]==tracker2[j] and tracker1[i]!='X': #find matched micelles
          dmf=comparison(steplist1[i],steplist2[j])
          totrf+=dmf
          mcount+=1

  return totrf,mcount

def infosplit(line):
  memline=line[96:]
  memline=memline.replace(",","")
  memsplit=memline.split()
  memsplit=[int(x) for x in memsplit]
  stsplit=line[6:18].split()
  site,time=int(stsplit[0]),float(stsplit[1])
  infoline=[site]+[time]+memsplit
  return infoline

#main fxn
def main():
  #Load input files
  morfile = open(sys.argv[1],'r')
  outfile = open(sys.argv[2],'w') #output file for kinetics summary
  #matrixfile = open(sys.argv[3],'w')
  outfile2= open(sys.argv[2]+'.rec','w')

  #set parameters if needed
  tstep=float(input("Time difference between snapshots (in ns)? ex) 0.02 (stop if irregular)\n"))
  nnmic=int(input("Most frequent number of micelles desired? ex) 30 or 8\n"))
  taumin=float(input("minimum correlation time you want (in ns)? ex) 0\n"))
  taumax=float(input("maximum correlation time you want (in ns)? ex) 100\n"))
  nmintau=int(taumin/tstep)
  nmaxtau=int(taumax/tstep)
  nsteptau=nmaxtau-nmintau

  start_time=timeit.default_timer()
  #read the clusterization result file, store all data in arrays
  #core array:3d list - [step][micelle][member]
  i,nstep=0,0
  biglist,steplist,nmiclist=[],[],[]
  for line in morfile:
    if line[0:4]=="step":#step initiating line
      ltsplit=line.split()
      nmcl1=int(ltsplit[3]) #number of mcls, temporary use
      nmiclist.append(nmcl1)
      i=0
      nstep+=1
    else: #collect
      i+=1 
      infolist=infosplit(line)
      miclist=infolist[2:]
      steplist.append(miclist)
      if i==nmcl1: #1 timestep ended
        biglist.append(steplist)
        steplist=[]
  elapsed=timeit.default_timer() - start_time
  print('finished reading cluster file, {} snapshots detected, time {:11.4f}'.format(nstep,elapsed))
  print('data size {} bytes'.format(sys.getsizeof(biglist)))

  #before calculate essential quantities, construct the tracker matrix
  tracker_m=[] #preparing first line of tracker matrix
  oldiline=[int(x) for x in range(1,nmiclist[0]+1)]
  tracker_m.append(oldiline)
  lastnorm=0 # #step which had normal #mic most recently
  for i in range(1,nstep):
    #matching adjacent steps: if #mic returns to normal value, try matching from the last step which
    #had normal number of micelles.
    #nmic1,nmic2=len(biglist[i-1]),len(biglist[i])
    #if nmic1!=nmic2 and nmic2==nnmic: 
    #  oldiline=tracker_m[lastnorm]
    #  matchline=tracker(biglist[lastnorm],biglist[i],oldiline)
    #else: matchline=tracker(biglist[i-1],biglist[i],oldiline)
    oldiline=tracker_m[lastnorm]
    matchline=tracker(biglist[lastnorm],biglist[i],oldiline)
    tracker_m.append(matchline)
    #oldiline=list(matchline)
    if len(matchline)==nnmic and 'X' not in matchline: lastnorm=i
  elapsed=timeit.default_timer() - start_time
  print('finished constructing tracker matrix, time {:11.4f}'.format(elapsed))
  #debug purpose:disclose the tracker matrix
  #for line in tracker_m:
  #  for x in line:
  #    matrixfile.write('{} '.format(x))
  #  matrixfile.write('\n')
  #matrixfile.close()

  #time correlation function loop
  #tau loop comes first to decide t loop. t loop is inner
  rf,count=numpy.zeros(nsteptau),numpy.zeros(nsteptau)
  rfrec,recount=numpy.zeros(nsteptau),numpy.zeros(nsteptau)
  taindex=0
  for j in range(nmintau+1,nmaxtau+1):
    for i in range(nstep-j):
      #collect molecule retention data of 1 type of tau
      rcorrinfo=relaxcorr_1step(biglist[i],biglist[i+j],tracker_m[i],tracker_m[i+j])
      rf[taindex]+=rcorrinfo[0]
      count[taindex]+=rcorrinfo[1]
      if nmiclist[i]==nnmic and nmiclist[i+j]==nnmic and 'X' not in (tracker_m[i]+tracker_m[i+j]):
        rfrec[taindex]+=rcorrinfo[0]
        recount[taindex]+=rcorrinfo[1]
    if j % 10 == 0:
      elapsed=timeit.default_timer() - start_time
      print('time correlating tau-step {}, time {:11.4f}'.format(j,elapsed))
    taindex+=1

  #output printing section
  outfile.write('Correlation_time(tau) Relaxation_corr(R(tau)) \n')
  for i in range(nsteptau):
    tau=taumin+(i+1)*tstep
    rf[i]/=count[i]
    rfrec[i]/=recount[i]
    outfile.write('{:9.4f} {:11.6f}\n'.format(tau,rf[i]))
    outfile2.write('{:9.4f} {:11.6f}\n'.format(tau,rfrec[i]))    

  morfile.close()
  outfile.close()
  outfile2.close()

if __name__ == "__main__": main()

