#!/usr/bin/env python
#program to read trajectory file and analyize ave box size 
#from npt simulation trajectory
#inputs: pdbfile, traj (dcd or xtc) file.

#import math
import sys
import numpy
import mdtraj as md

#def find_nearest_index(array, value):
    #array = np.asarray(array)
#    idx = (numpy.abs(array - value)).argmin()
#    return idx

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  #outfile = sys.argv[3]
  #boxdemand=float(input("How much is the box size demanded(in nm)? ex) 4.34 \n"))
  teq=int(input("How many initial frames do you want to cut as equilibration? ex) 5000 \n"))
  tend=int(input("Until which frame(put -1 for limitless)? ex) 10000 \n"))

  #input 1 : load surf traj. (big file)
  traj=md.load(trjfile,top=topfile)
  if tend==-1:
    traj=traj[teq:]
  else:
    traj=traj[teq:tend]
  nstep=traj.n_frames
  boxinfo=traj.unitcell_lengths 
  #traj.save_xtc(outfile)

  numpy.savetxt('boxdim'+trjfile[:-4]+'.dat',boxinfo)
  avedim=numpy.mean(boxinfo,axis=0)
  print(avedim)
  #idx=find_nearest_index(boxinfo,boxdemand)
  #snapshot=traj[idx]
  #topology=snapshot.topology
  #drude particle exclusion
  #snapshot=snapshot.atom_slice([atom.index for atom in topology.atoms if ('D' not in atom.name) ])
  #snapshot.save_pdb(outfile)

if __name__ == "__main__": main()

