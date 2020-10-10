#!/usr/bin/env python
#instant program to fix issue 

import mdtraj as md
import numpy as np

ndim=2000
instr="Urea_Urea_"
desti="image_processed/"

#def swap_coordinates(trajectory,i,j):
#top=t.topology
#coordinates=t.xyz
#coordinates[0][j],coordinates[0][i]=coordinates[0][i],coordinates[0][j]
#swap_coordinates(t,0,1)

tottraj=md.load("Urea_Urea_0.pdb",top="Urea_Urea_0.pdb")

for i in range(ndim):
  pdbin=instr+str(i)+".pdb"
  pdbout=desti+pdbin
  xyzout=desti+instr+str(i)+".xyz"
  traj=md.load(pdbin,top=pdbin)

  bondsarray=np.array([ [0,10],[10,8],[8,9],[8,11],[9,14],[9,15],[11,12],[11,13]  ])

  traj.make_molecules_whole(sorted_bonds=bondsarray)
  traj=traj.center_coordinates()

  traj[0].save_pdb(pdbout)
  traj[0].save_xyz(xyzout)
  if i==0:
    tottraj=traj
  else:
    tottraj=tottraj.join(traj)

tottraj.save_pdb(desti+"tottraj_urea_urea.pdb")
tottraj.save_xyz(desti+"tottraj_urea_urea.xyz")

