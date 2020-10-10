#!/home/kjeong23/softwares/bin/python3.4
# micelle kinetics analysis program
# input file : merged lattice point assignment summary file for entire time period
#v02: site-specific kinetic info, micelle exchange portion, exchange rate, gain&lose
#v03: retention time histogram, raw data record for (m-index,psite#,site#,ret.time)
#v03 renovation: also counted diffusion out to solvent media.

import math
import sys
import numpy

def comparison(a1,a2): #getting retention fraction of micelle during change a1->a2
  coin=0
  for x in a1:
    if x in a2:
      coin+=1
  retf=coin/len(a1)
  return retf

#main fxn
def main():
 
  #Load input files
  morfile = open(sys.argv[1],'r')
  outfile = open(sys.argv[2],'w') #output file for kinetics summary

  nsurf=int(input("Total # of surfactant molecules? ex) 973 \n"))
  nlat=int(input("Total # of lattice points in the system? ex) 30 \n"))
  begt=float(input("Initial time of the whole trajectory in ns? ex) 300.00 \n"))
  fint=float(input("Final time of the whole trajectory in ns? ex) 499.96 \n"))
  #stept=float(input("Snapshot recording timestep in ns? ex) 0.040 \n"))
  bs=float(input("Bin size of hop ret.time histogram? ex) 1\n"))
  nbin=int( math.ceil((fint-begt)/bs) ) # of bins
  psit=numpy.zeros(nbin) #psi(t). probability distribution of hopping retention time
  rawret=numpy.empty((0,4),float) # raw data storage for retention time

  latn=[0]*nsurf #array of lattice point index of each surf
  platn=[0]*nsurf #array of past lattice point indices
  stackrett=numpy.zeros(nsurf) #stacked time array for retention of each surf
  initt=[-1]*nsurf #array of initial time of retention of each surf
  lastt=[-1]*nsurf #array of final time of retention of each surf
  rett=[0.0]*nsurf #total retention time regardless of micelle#
  homem=[1]*nsurf #number of home micelles in the interval, for each surf
  imem,fmem=[[]]*nlat,[[]]*nlat #array of initial,final members of each lattice point
  gain,lose=[0]*nlat,[0]*nlat #counting surfactant gaining,losing of each lattice point

  ovhop,i,ii=0,0,0
  #start loop of reading. renov: impose a 'cycle' for each timestep
  for line in morfile:
    if len(line)>=92: #for the line actually telling members
      i+=1
      #ii+=1
      memline=line[92:]
      memline=memline.replace(",","")
      memsplit=memline.split() #1d list of member surf indices
      #if ii<=30:
      #  print(memsplit)
      ltsplit=line[6:17].split() #lattice point index and time
      site,time=int(ltsplit[0]),float(ltsplit[1])
      if imem[site-1]==[]: #first step (data not registered at all)
        imem[site-1]=memsplit #register initial members
      if time==fint: #final step
        fmem[site-1]=memsplit #register final members
      for x in memsplit: #examination for each member surfactant
        y=int(x)
        if time==begt: #initial step
          #print('y',y)
          #print('site',site)
          latn[y],platn[y]=site,site #unregistered(in solvent media): #0
          ptime=time #past step's time
        else: #not the 1st step
          latn[y]=site
    if i>=nlat: #end of one snapshot
      i=0
      for j in range(nsurf):
        if latn[j]==platn[j]: #lattp maintained
          if latn[j]!=0:
            stackrett[j]+=(time-ptime) #stack the retention time
        else: #lattp changed
          if platn[j]!=0: #not from solvent region
            #store stacked result into permanent record and histogram
            onestack=numpy.array( [ j , platn[j] ,latn[j], stackrett[j] ] )
            rawret=numpy.vstack((rawret,onestack))
            psit[ int( math.ceil(stackrett[j]/bs) )-1]+=1
            ovhop+=1 
            stackrett[j]=0 #initialize
            #lose[platn[j]-1]+=1
            #gain[site-1]+=1
            #homem[j]+=1
      #at the end of whole trajectory: conclude stackable data
      if time==fint:
        for j in range(nsurf):
          if latn[j]!=0 and stackrett[j]!=0: #don't overcount 'newly changed' at last step
            onestack=numpy.array( [ j , platn[j] ,latn[j], stackrett[j] ] )
            rawret=numpy.vstack((rawret,onestack))
            psit[ int( math.ceil(stackrett[j]/bs) )-1]+=1
            ovhop+=1
      #at the end of each snapshot: initialize 'current latn array' to 0(solvent)
      #if not in solvent, this will be overwritten anyhow
      for j in range(nsurf):
        platn[j]=latn[j]
        latn[j]=0
      ptime=time #update the past step's time

  #after end of loop, final process 1. hopping statistics
#  ovave, ovhop=0.0, 0
#  for y in range(nsurf):
#    rett[y]+=(lastt[y]-initt[y]) #counting the retention in last home micelle
#    outfile.write('Surf {:3} Ave ret.time {:7.3f} ns, #hopping {} \n'.format(y,rett[y]/homem[y],homem[y]-1))
#    ovave+=rett[y]/homem[y]
#    ovhop+=homem[y]-1
#  ovave/=nsurf
#  avhop=ovhop/nsurf
#  outfile.write('Overall average retention time per surfactant : {} ns \n'.format(ovave))
#  outfile.write('Tot# of hopping {} times during whole traj, average {:7.3f} times per molecule \n'.format(ovhop,avhop))

  #micelle retention fraction
#  for x in range(nlat):
#    retf=comparison(imem[x],fmem[x])
#    outfile.write('lattp# {:2} mem retention frtn {:6.3f} mem# {:2} -> {:2} receive {:2} give {:2} \n'.format(x+1,retf,len(imem[x]),len(fmem[x]),gain[x],lose[x]))

  #histogram section
  psit=(psit/ovhop)/bs #normalize
  outfile.write('Retention time histogram: tau(ps) psi(tau)\n')
  for x in range(nbin):
    outfile.write('{:8.4f} {:8.4f}\n'.format((x+1)*bs,psit[x]))

  #raw data section
  outfile.write('Raw data: Surf#, site#, ret.time. note: total hopping {:6}\n'.format(ovhop))
  for line in rawret:
    outfile.write('{:4} {:3} {:3} {:8.4f}\n'.format(int(line[0]),int(line[1]),int(line[2]),line[3]))

  morfile.close()
  outfile.close()

if __name__ == "__main__": main()

