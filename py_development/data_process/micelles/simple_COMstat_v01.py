#!/home/kjeong23/softwares/bin/python3.4
# micelle COM distance tracker. 
# input file : entire COM file, migration event file (should get rid of explanations) 
# output file : distance distribution of contact vs non-contact
# be careful that micelle number in input files starts from 1. So for array index, I have to subtract 1

import sys
import numpy
import timeit

def arrayread_com_1step(list_2d,nmic):
  com_array=numpy.full((nmic,3),-1000) #dummy value to notate unassigned micelle
  for row in list_2d:
    rowsplit=row.split()
    mic=int(rowsplit[1])-1
    com_array[mic]=[float(rowsplit[4]),float(rowsplit[5]),float(rowsplit[6])]
  return com_array

def pbcdist_3d(r1,r2,box): #distance regarding pbc
  r=r2-r1
  for i in range(3):
    if r[i]>(box[i]/2.0):
      r[i]-=box[i]
    elif r[i]<(-box[i]/2.0):
      r[i]+=box[i]
  dist=numpy.sqrt(numpy.dot(r,r))
  return dist

def main():
  #Load input files
  comfile = open(sys.argv[1],'r')
  migfile = open(sys.argv[2],'r')
  outfile = open(sys.argv[3],'w') #output file name.(COM dist. distrib. for exch,non-exch)

  #nmic=int(input("Most frequent number of micelles desired? ex) 30 or 8\n"))
  nmic=30
  #boxinfo=input("3d periodic Box dimensions? ex) 13.508 13.508 7.104")
  #box=numpy.array([float(x) for x in boxinfo.split()])
  box=numpy.array([13.508,13.508,7.104])
  #tstep=float("timestep between recorded snapshots (in ns)? ex) 0.02")
  tstep=0.02
  #nstep=int("Total number of snapshots? ex) 50000")
  nstep=50000
  micpair=input("Which pair of micelle sites do you want to inspect? ex) 2 20\n")
  mic12=micpair.split()
  mic1,mic2=int(mic12[0]),int(mic12[1])
  distinfo = input("Min,Max and Grid size for Prob.dist of COM distance? ex) 0 5 0.02\n")
  dset =distinfo.split()
  dmin,dmax,dbin=float(dset[0]),float(dset[1]),float(dset[2])
  nbin=int((dmax-dmin)/dbin)+1
  start_time=timeit.default_timer()

  sindex=0
  com=numpy.empty((nstep,nmic,3))
  state=numpy.ones(nstep) #distinguish whether a snapshot is exchanging surfactant or not
  #d=numpy.empty(nbin)
  #for i in range(nbin):
  #  d[i]=dmin+dbin*i
  pdd1,pdd2=numpy.zeros(nbin),numpy.zeros(nbin) #non-exch COMdist distrib. and exch.
  #read COM information, compose matrix
  while True:
    comline=comfile.readline()
    if comline=='': #EOF
      break
    if comline[0:4]=="step":#step initiating line
      ltsplit1=comline.split()
      nmcl1=int(ltsplit1[3])
      paragraph=[]
      for i in range(nmcl1):
        comstr=comfile.readline()
        paragraph.append(comstr)
      com[sindex]=arrayread_com_1step(paragraph,nmic)
    sindex+=1

  elapsed=timeit.default_timer() - start_time
  print('COM coordinate reading complete. Time {}'.format(elapsed))

  #read exchange event file, distinguish exchange step and non-exchange step
  #at this stage, not specifically distinguishing acc or don.
  #complete the state file. non-exchange:1, exchange:2
  while True:
    migline=migfile.readline()
    if migline=='':
      break
    migsplit=migline.split()
    accm,donm=float(migsplit[1]),float(migsplit[2])
    accm,donm=int(round(accm)),int(round(donm))
    if [accm,donm]==[mic1,mic2] or [accm,donm]==[mic2,mic1]: #detected corresponding exch
      step1,step2=float(migsplit[3])/tstep,float(migsplit[4])/tstep
      step1,step2=int(step1),int(step2)
      for i in range(step1,step2+1):
        state[i]=2

  #COM distance
  pymic1,pymic2=mic1-1,mic2-1
  for i in range(nstep):
    r1,r2=com[i][pymic1],com[i][pymic2]
    if r1[0]!=-1000 and r2[0]!=-1000:
      ccd=pbcdist_3d(r1,r2,box)
      index=int((ccd-dmin+dbin/2.0)/dbin)
      if index<nbin:
        if state[i]==1:
          pdd1[index]+=1
        elif state[i]==2:
          pdd2[index]+=1

  #normalize and writing section
  print(sum(pdd1),sum(pdd2))
  pdd1=pdd1/(numpy.sum(pdd1)*dbin)
  pdd2=pdd2/(numpy.sum(pdd2)*dbin)
  for i in range(nbin):
    outfile.write('{:8.3f} {:8.3f} {:8.3f}\n'.format(dmin+dbin*i,pdd1[i],pdd2[i]))

  comfile.close()
  migfile.close()
  outfile.close()
   
if __name__ == "__main__": main()
