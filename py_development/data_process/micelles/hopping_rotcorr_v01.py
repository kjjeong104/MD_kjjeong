#!/home/kjeong23/softwares/bin/python3.4
# surfactant hopping rotational correlation function calculator
# algorithm : get coord.(big file)& user input -> starts loop of migevent reading.
# -> for each relevant migration event, cut trajectory of 1 surf by x ns from departure time.
# -> calc rot.corr fxn and collect. calc average by migration events later.
# input : grofile trjfile migfile, output : rot corr fxn file.
# (User input): (range of tau), trajectory trimming length for event, time-block information.

import math
import sys
import numpy
import timeit
import mdtraj as md

tstep=0.02

def mig_interpret(migline,tb1,tb2,trimt): #read migline, interpret if this hopping event is relevant
  #read surf index and initial time(important)
  split=migline.split()
  mindex,t1=split[0],float(split[3]) #mindex : 'string' form
  t2=t1+trimt
  decision=False  
  if tb1<=t1 and t2<=tb2: decision=True
  hop_info=[decision,mindex,t1-tb1,t2-tb1]
  return hop_info

def rot_acf_lmol(trimtraj,i1,i2): #calculator for rotational acf of a linear molecule
  #make sure that molecules in the trajectory are whole. trimtraj:1molecule traj during an interval
  #i1, i2 : atom indices for two points to indicate orientation.
  #rot_acf : maybe numpy array.
  nstep,crd=trimtraj.n_frames,trimtraj.xyz
  nsteptau=int(nstep/2)+1
  rot_acf=numpy.zeros(nsteptau)
  #loop of autocorrelation
  for j in range(nsteptau):
    count=0
    for i in range(nstep-j):
      v0,vt=crd[i][i2]-crd[i][i1],crd[i+j][i2]-crd[i+j][i1] #pick vectors to self-correlate
      v0,vt=v0/numpy.linalg.norm(v0),vt/numpy.linalg.norm(vt) #orientations: set to unit vector
      rot_acf[j]+=numpy.dot(v0,vt)
      count+=1
    rot_acf[j]/=count #normalize 

  return rot_acf

def main():
  #input file (trajfile: binary processing)
  grofile = sys.argv[1]
  trjfile = sys.argv[2]
  migfile = open(sys.argv[3],'r')
  outfile = open(sys.argv[4],'w')

  i1,i2=0,10 #internal atom indices in a surfactant molecule which directs tail end, headgroup P
  trimt=float(input("migration trajectory trimming time you want (in ns)? ex) 2\n"))
#  taumin=float(input("minimum correlation time you want (in ns)? ex) 0\n"))
#  taumax=float(input("maximum correlation time you want (in ns)? ex) 100\n"))
#  nmintau=int(taumin/tstep)
#  nmaxtau=int(taumax/tstep)
#  nsteptau=nmaxtau-nmintau
  nsteptau = int(trimt/tstep/2)+1
  print(nsteptau)
## when number of frames are not big, then we don't have to divide tau interval to do HT computing.
  tblock=input("What is the time interval of trajectory you're loading(in ns)? ex) 200 400 \n")
  tb12=tblock.split()
  tb1,tb2=float(tb12[0]),float(tb12[1])
  start_time=timeit.default_timer()

  #input 1 : load surfactant trajectory.(big file)
  traj=md.load(trjfile,top=grofile)
  topology=traj.topology
  #print(traj)
  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))

  #input 2 & loop : load migevents, and for each migevent, cut traj and calc rotacf on the fly.
  lindex,cfcount=0,0
  tot_rtcf=numpy.zeros(nsteptau)
  while True:
    migline=migfile.readline()
    if migline=='':#EOF
      break
    info=mig_interpret(migline,tb1,tb2,trimt)
    if info[0]==True: #decided to include that hopping event in calculation
      nt1,nt2=info[2]/tstep,info[3]/tstep
      nt1,nt2=int(round(nt1)),int(round(nt2))
      trimtraj=traj[nt1:nt2+1]
      trimtraj=trimtraj.atom_slice(topology.select('resid '+info[1]))
      onemol_rtcf=rot_acf_lmol(trimtraj,i1,i2)
      tot_rtcf+=onemol_rtcf
      cfcount+=1
    if lindex%1000==0:
      elapsed=timeit.default_timer() - start_time
      print('examining migration event # {}, time {:11.4f}'.format(lindex,elapsed))
    lindex+=1

  #output writing section
  tot_rtcf/=cfcount#normalize
  for i in range(tot_rtcf.size):
    tau=i*tstep
    outfile.write('{:9.4f} {:11.6f}\n'.format(tau,tot_rtcf[i]))

  migfile.close()
  outfile.close()

if __name__ == "__main__": main()
