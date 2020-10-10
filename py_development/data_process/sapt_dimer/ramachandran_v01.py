##python script for peptide ramachandran plot
# molecule should be in unwrapped trajectory.
import sys
import math
import mdtraj as md
#from __future__ import print_function
#from simtk.openmm.app import *
#from simtk.openmm import *
#from simtk.unit import *
#from sys import stdout
#from time import gmtime, strftime
#from datetime import datetime
import numpy as np

radtodeg=180.0/math.pi
pdbin=sys.argv[1]
trjin=sys.argv[2]
outfile=open(sys.argv[3],'w')
phiin=input("pdbfile atom indices for phi dihedral angle. ex) 8 7 9 11 \n")
psiin=input("pdbfile atom indices for psi dihedral angle. ex) 16 15 9 11 \n")
#construct atom index array (need to shift by -1, in order to  )
phiind=np.array([[int(x)-1 for x in phiin.split()]])
psiind=np.array([[int(x)-1 for x in psiin.split()]])

traj=md.load(trjin,top=pdbin)
topology=traj.topology
nstep=traj.n_frames
phi=md.compute_dihedrals(traj,phiind,periodic=True)
psi=md.compute_dihedrals(traj,psiind,periodic=True)

for i in range(nstep):
  outfile.write("{:8.3f} {:8.3f}\n".format(phi[i][0]*radtodeg,psi[i][0]*radtodeg))

outfile.close()

