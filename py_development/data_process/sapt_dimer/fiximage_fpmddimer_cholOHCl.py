#!/usr/bin/env python
#instant program to fix issue 

import mdtraj as md
import numpy as np

ndim=1000
instr="Chol_Cl_"
desti="image_processed/"

#def swap_coordinates(trajectory,i,j):
#top=t.topology
#coordinates=t.xyz
#coordinates[0][j],coordinates[0][i]=coordinates[0][i],coordinates[0][j]
#swap_coordinates(t,0,1)

tottraj=md.load("Chol_Cl_0.pdb",top="Chol_Cl_0.pdb")

#for i in range(710,711):
for i in range(ndim):
  pdbin=instr+str(i)+".pdb"
  pdbout=desti+pdbin
  xyzout=desti+instr+str(i)+".xyz"
  oldtraj=md.load(pdbin,top=pdbin)
  oldtopol=oldtraj.topology
  #Cl swap recovery part
  if 'Cl' in str(oldtopol.atom(0)): #first atom is Cl -> chol cl swapped
  #then among the 22 atoms, index 0(Cl) -> 21. chol atoms 1~21 -> 0~20.
    print(i)
    traj1=oldtraj.atom_slice([0]) #Cl
    #print(traj1)
    traj2=oldtraj.atom_slice(range(1,22)) #rest
    #print(traj2)
    traj=traj2.stack(traj1)
    print(traj)
    #print(traj)
    #topology=traj.topology
    #traj.xyz[0][21]=oldtraj.xyz[0][0] #move Cl to the end
    #traj.xyz[0][0:21],traj.xyz[0][21]=oldtraj.xyz[0][1:22],oldtraj.xyz[0][0] #move Cl to the end
    #topology.atom(range(0,21)),topology.atom(21)=oldtopol.atom(range(1,22)),oldtopol.atom(0)
    #for k in range(0,21):
    #  traj.xyz[0][k]=oldtraj.xyz[0][k+1]
  else:
    traj1,traj2=oldtraj.atom_slice(range(0,21)),oldtraj.atom_slice([21])
    traj=traj1.stack(traj2)
  #unwrapping part
  bondsarray=np.array([ [20,21] ])

  traj.make_molecules_whole(sorted_bonds=bondsarray)
  traj=traj.center_coordinates()

  traj[0].save_pdb(pdbout)
  traj[0].save_xyz(xyzout)
  if i==0:
    tottraj=traj
  else:
    tottraj=tottraj.join(traj)

tottraj.save_pdb(desti+"tottraj_cholOHCl.pdb")
tottraj.save_xyz(desti+"tottraj_cholOHCl.xyz")

