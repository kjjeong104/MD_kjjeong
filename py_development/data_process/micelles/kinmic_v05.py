#!/home/kjeong23/softwares/bin/python3.4
# micelle kinetics analysis program
# input file : merged/reassigned lattice point assignment summary file for entire time period
# doesn't need tracker file now, but the assignment should be already done properly.
#v02: site-specific kinetic info, micelle exchange portion, exchange rate, gain&lose
#v03: retention time histogram, raw data record for (m-index,psite#,site#,ret.time)
#v03 renovation: also counted diffusion out to solvent media.
#v04:now should be able to count even if clustering data has dislocation/fission.
#also calculate "time dependent relaxation function" of member retention
#policy for fission/fusion : just "uncounting" the broken/coalesced micelles
#v05: count exchange events and write log. (reviving v02,v03,m-index,surf#,psite#,nsite#,arrival time)

import math
import sys
import numpy
import timeit

def comparison(a1,a2): #getting retention fraction of micelle during change a1->a2
  coin=0
  for x in a1:
    if x in a2:
      coin+=1
  retf=coin/len(a1)
  return retf

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
  outfile = open(sys.argv[2],'w') #output file for event history

  nnmic=int(input("Most frequent number of micelles desired? ex) 30 or 8\n"))
  nstep=int(input("Total number of steps? ex) 50000\n")) 

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

  #event counting


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

