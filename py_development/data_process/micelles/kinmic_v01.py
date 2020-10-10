#!/home/kjeong23/softwares/bin/python3.4
# micelle kinetics analysis program
# input file : merged lattice point assignment summary file for entire time period
#v01: average retention time for a single surfactant

import math
import sys

#main fxn
def main():
 
  #Load input files
  morfile = open(sys.argv[1],'r')
  outfile = open(sys.argv[2],'w') #output file for kinetics summary

  nsurf=int(input("Total # of surfactant molecules? ex) 973 \n"))

  latn=[0]*nsurf #array of lattice point index of each surf
  initt=[-1]*nsurf #array of initial time of retention of each surf
  lastt=[-1]*nsurf #array of final time of retention of each surf
  rett=[0.0]*nsurf #total retention time regardless of micelle#
  homem=[1]*nsurf #number of home micelles in the interval, for each surf

  #start loop of reading.
  for line in morfile:
    if len(line)>=56: #for the line actually telling members
      memline=line[58:]
      memline=memline.replace(",","")
      memsplit=memline.split() #1d list of member surf indices
      ltsplit=line[6:17].split() #lattice point index and time
      ltsplit[0],ltsplit[1]=int(ltsplit[0]),float(ltsplit[1])
      for x in memsplit: #examination for each member surfactant
        y=int(x)
        if latn[y]==0: #initial step
          latn[y],initt[y],lastt[y]=ltsplit[0],ltsplit[1],ltsplit[1]
        elif latn[y]==ltsplit[0]: #lattp maintained 
          lastt[y]=ltsplit[1]
        elif latn[y]!=ltsplit[0]: #lattp changed
          rett[y]+=(lastt[y]-initt[y])
          homem[y]+=1
          initt[y],lastt[y]=ltsplit[1],ltsplit[1]
          latn[y]=ltsplit[0]

  #after end of loop, final process
  ovave=0.0
  ovhop=0
  for y in range(nsurf):
    rett[y]+=(lastt[y]-initt[y]) #counting the retention in last home micelle
    outfile.write('Surf {} Ave ret.time {:7.3f} ns, #hopping {} \n'.format(y,rett[y]/homem[y],homem[y]-1))
    ovave+=rett[y]/homem[y]
    ovhop+=homem[y]-1
  ovave/=nsurf
  avhop=ovhop/nsurf
  outfile.write('Overall average retention time per surfactant : {} ns \n'.format(ovave))
  outfile.write('Tot# of hopping {} times during whole traj, average {:7.3f} times per molecule \n'.format(ovhop,avhop))

  morfile.close()
  outfile.close()

if __name__ == "__main__": main()

