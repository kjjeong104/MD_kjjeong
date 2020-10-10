#!/home/kjeong23/softwares/bin/python3.4
# program to calculate rotational correlation function of an arbitrary molecule (SAPT or GAFF)
# algorithm : read dcd file to get coord -> cut traj of 1 molecule -> calc rot.corr fxn
# -> collect, calc average
#inputs: pdb file (topology topfile), traj dcd file. output: rot corr fxn file.

import math
import sys
import numpy
import timeit
import mdtraj as md

def rot_acf_lmol(trimtraj): #calculator for rotational acf of a linear molecule
  #initial-time averaged, and initial-time not averaged(only at original time)
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
      v0,vt=crd[i][1]-crd[i][0],crd[i+j][1]-crd[i+j][0] #pick vectors to self-correlate
      v0,vt=v0/numpy.linalg.norm(v0),vt/numpy.linalg.norm(vt) #orientations: set to unit vector
      rot_acf[j]+=numpy.dot(v0,vt)
      count+=1
    rot_acf[j]/=count #normalize 
  return rot_acf

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  outfile = open(sys.argv[3],'w')

  aname=input("What is the atom pair name to define orientation? ex) O2a C2a or Nm1 Nm2 \n")
  an12=aname.split()
  an1,an2=str(an12[0]),str(an12[1])
  #tblock=input("What is the time interval of trajectory you're loading(in ns)? ex) 200 400 \n")
  #tb12=tblock.split()
  #tb1,tb2=float(tb12[0]),float(tb12[1])
  #ntb1,ntb2=tb1/tstep,tb2/tstep
  #ntb1,ntb2=int(round(ntb1)),int(round(ntb2))
  mi12=input("What is the residue number index interval of rotacf calculation? ex) 0 299 \n")
  mi12=mi12.split()
  mi1,mi2=int(mi12[0]),int(mi12[1])
  tstep=float(input("timestep between snapshots in ps? ex) 0.01\n"))
  tskip=int(input("Once in how many frames do you want to take? ex) 1 \n"))

  start_time=timeit.default_timer()

  #input 1 : load surf traj. (big file)
  traj=md.load(trjfile,top=topfile)
  traj=traj[::tskip]
  #tstep=traj.timestep
  topology=traj.topology
  nstep=traj.n_frames
  nsteptau=int(nstep/2) + 1
  nmon=topology.n_residues
  #traj=traj.atom_slice(topology.select('name =='+str(aname)))
  if nmon <= mi2: #residue index exceeded
    mi2=nmon-1 #autofix surfactant index range
  #boxinfo=numpy.empty((0,3)) #N*3 array notates box dimension vector for all snapshots
  #for i in range(nstep):
  #  boxone=numpy.array(traj.unitcell_lengths[i])
  #  boxinfo=numpy.vstack((boxinfo,boxone))

  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))
  print(nstep,nmon,nsteptau)

  #prepare bins for rotacf
  tot_rtcf=numpy.zeros(nsteptau)
  cfcount=0

  #loop of trajectory slicing and rotacf calculation
  for i in range(mi1,mi2+1):
    trimtraj=traj.atom_slice(topology.select('resid '+str(i)))
    topone=trimtraj.topology
    trimtraj=trimtraj.atom_slice(topone.select('(name == '+str(an1)+') or (name == '+str(an2)+')'))
    #trimtraj=trimtraj.atom_slice(topology.select('name =='+str(aname)))
    onemol_rtcf=rot_acf_lmol(trimtraj)
    tot_rtcf+=onemol_rtcf
    cfcount+=1

    elapsed=timeit.default_timer() - start_time
    print('molecule# {} rotacf calculated. time{:11.4f}'.format(i,elapsed))

  #output writing section
  tot_rtcf/=cfcount#normalize
  for i in range(nsteptau):
    tau=float(i)*tstep*tskip
    outfile.write('{:9.4f} {:11.6f}\n'.format(tau,tot_rtcf[i]))

  outfile.close()

if __name__ == "__main__": main()

