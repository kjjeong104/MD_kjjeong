#!/home/kjeong23/softwares/bin/python3.4
# program to calculate MSD of an arbitrary molecule (SAPT or GAFF)
# algorithm : read dcd file to get coord -> unwrapping coordinate and calc MSD.
# for simplicity purpose, pick C atom of urea.
# imported mdtraj to process dcd file. dcd traj is truncated with 1 atomtype.
# cut traj. -> record r0 (don't calc COM) 
# in every step of COM SD calc, pile up in the bins, adding up "counts"(w.r.t. t)
# to compare statistical weight of all dr^2(t). -> calc MSD(t).
# -> display MSD(t), count(t).
#inputs: pdbfile, traj dcd file
#output: file of {t, MSD(t), count(t)}
#v011 : able to cut some interval initially for equilibration run

import math
import sys
import numpy
import timeit
import mdtraj as md

#tstep=1.000 #used ps time unit. 1 snapshot : 1ps.

def pbcdr(ri,rf,box): #displacement vector with considering pbc
  dr=rf-ri
  for i in range(3):
    if dr[i]>(box[i]/2.0):
      dr[i]-=box[i]
    elif dr[i]<(-box[i]/2.0):
      dr[i]+=box[i]
  return dr

def SD_calculator(trimtraj,nstep,boxinfo): #square displacement calculator. 
  #assumes that trajectory has 1 particle. however, box can change (in case of NPT)
  sdbin1,count1=numpy.zeros(nstep),numpy.zeros(nstep)
  trimstep=trimtraj.n_frames
  crd=trimtraj.xyz
  #box=numpy.array(trimtraj.unitcell_lengths[0])
  #for ntau in range(trimstep): #use tau-step loop first
  #  for i in range(trimstep-ntau): #initial snapshot
  for i in range(trimstep-1):
    vec=numpy.zeros(3)
    for j in range(i+1,trimstep):
      dr=pbcdr(crd[j][0],crd[j-1][0],boxinfo[j])
      vec+=dr
      rmag2=numpy.dot(vec,vec)
      sdbin1[j-i]+=rmag2
      count1[j-i]+=1
  sdbin1[0],count1[0]=0.0,1
  return sdbin1,count1

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  outfile = open(sys.argv[3],'w')

  aname=input("What is the atom name to track MSD? ex) C2a \n")
  #tblock=input("What is the time interval of trajectory you're loading(in ns)? ex) 200 400 \n")
  #tb12=tblock.split()
  #tb1,tb2=float(tb12[0]),float(tb12[1])
  #ntb1,ntb2=tb1/tstep,tb2/tstep
  #ntb1,ntb2=int(round(ntb1)),int(round(ntb2))
  mi12=input("What is the residue number index interval of MSD calculation? ex) 0 299 \n")
  mi12=mi12.split()
  mi1,mi2=int(mi12[0]),int(mi12[1])
  tskip=int(input("Once in how many frames do you want to take? ex) 10 \n"))
  teq=int(input("How many initial frames do you want to cut as equilibration? ex) 5000 \n"))

  start_time=timeit.default_timer()

  #input 1 : load surf traj. (big file)
  traj=md.load(trjfile,top=topfile)
  traj=traj[teq::tskip]
  tstep=float(traj.timestep)
  topology=traj.topology
  nstep=traj.n_frames
  nmon=topology.n_residues
  #traj=traj.atom_slice(topology.select('name =='+str(aname)))
  if nmon <= mi2: #residue index exceeded
    mi2=nmon-1 #autofix surfactant index range
  boxinfo=numpy.empty((0,3)) #N*3 array notates box dimension vector for all snapshots
  for i in range(nstep):
    boxone=numpy.array(traj.unitcell_lengths[i])
    boxinfo=numpy.vstack((boxinfo,boxone))

  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))
  print(nstep,nmon)

  #prepare bins for MSD
  sdbin,count=numpy.zeros(nstep),numpy.zeros(nstep) #bin & stat weight of dr^2(t) ave

  #loop of trajectory slicing and MSD calculation
  for i in range(mi1,mi2+1):
    trimtraj=traj.atom_slice(topology.select('resid '+str(i)+' and name == '+str(aname)))
    #trimtraj=traj.atom_slice(topology.select('name == '+str(aname)))
    #newtop=trimtraj.topology
    #print(trimtraj,newtop)
    #for res in newtop.residues:
    #  print(res.index)
    #print(newtop.select('resid '+str(i)))
    #trimtraj=trimtraj.atom_slice(newtop.select('resid '+str(i)))
    #trimtraj=trimtraj.atom_slice(topology.select('name =='+str(aname)))
    sdbin1,count1=SD_calculator(trimtraj,nstep,boxinfo)
    sdbin+=sdbin1
    count+=count1

    elapsed=timeit.default_timer() - start_time
    print('molecule# {} MSD calculated. time{:11.4f}'.format(i,elapsed))

  #printing section
  sdbin=numpy.divide(sdbin,count,out=numpy.zeros_like(sdbin), where=count!=0)
  #outfile.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(0,0,0))
  for i in range(nstep):
    outfile.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(tstep*float(i),sdbin[i],count[i]))

  outfile.close()

if __name__ == "__main__": main()

