#!/usr/bin/env python
#instant program to fix issue 

import mdtraj as md
import numpy as np

ndim=1000
instr="Chol_Urea_"
desti="image_processed/"

tottraj=md.load("Chol_Urea_0.pdb",top="Chol_Urea_0.pdb")

for i in range(ndim):
  pdbin=instr+str(i)+".pdb"
  pdbout=desti+pdbin
  xyzout=desti+instr+str(i)+".xyz"
  traj=md.load(pdbin,top=pdbin)

  bondsarray=np.array([ [20,23],[23,21],[21,22],[21,24],[22,27],[22,28],[24,25],[24,26]  ])

  traj.make_molecules_whole(sorted_bonds=bondsarray)
  traj=traj.center_coordinates()

  traj[0].save_pdb(pdbout)
  traj[0].save_xyz(xyzout)
  if i==0:
    tottraj=traj
  else:
    tottraj=tottraj.join(traj)

tottraj.save_pdb(desti+"tottraj_cholOHuO.pdb")
tottraj.save_xyz(desti+"tottraj_cholOHuO.xyz")

