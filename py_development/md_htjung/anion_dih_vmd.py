#!/usr/bin/env python2
# ver 0.1 - coding python by Hyuntae Jung on 5/13/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='visualize cations near cisoid TFSI using B-factor in PDB')
## args
parser.add_argument('-i', '--input', default='md.dcd', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='md.pdb', nargs='?', 
	help='pdb file for structure')
parser.add_argument('-d1', '--data1', default='cation.raw.dih.npy', nargs='?', 
	help='first dihedrals for each group - usually cation')
parser.add_argument('-d2', '--data2', default='anion.raw.dih.npy', nargs='?', 
	help='second dihedrals for each group - usually anion')
parser.add_argument('-b', '--begin', default=0, nargs='?', type=int,
	help='beginning iframe to visualize with cisoid conformation of TFSI')
parser.add_argument('-e', '--end', default=-1, nargs='?', type=int,
	help='end iframe to visualize with cisoid conformation of TFSI (-1 for last frame)')
parser.add_argument('-step', '--step', default=1, nargs='?', type=int,
	help='time interval between frames you want to save')
parser.add_argument('-o', '--output', default='cis_tfsi_vmd', nargs='?', 
	help='output prefix')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

args.outpdb = args.output + '.pdb'

## import modules
import mdtraj as md
import MDAnalysis as mda
import numpy as np
import math

## load dihedral angles
dih_cat = np.load(args.data1)
dih_cat = dih_cat*180.0/math.pi/10.0 # reduced by 10 for writing pdb file
dih_ani = np.load(args.data2)
dih_ani = dih_ani*180.0/math.pi/10.0 # reduced by 10 for writing pdb file
dih = np.append(dih_cat,dih_ani,axis=1)
print(" assume cation is placed before anion in topology file")

# read 
u = mda.Universe(args.structure, args.input)
real_atoms = u.select_atoms("all and not name D*")
n_residues = len(real_atoms.residues.indices)
if dih.shape[1] != n_residues:
	raise ValueError("total number of residue is not matched. {} residues != {} dihedrals".format(n_residues,dih.shape[1]))
n_frames = len(u.trajectory)

if args.end == -1:
	args.end = n_frames
elif args.end < args.begin:
	raise ValueError("{} end frame is less than {} begin frame?".format(args.end,args.begin))
elif args.end > n_frames:
	raise ValueError("{} end frame is greater than {} n_frames?".format(args.end,n_frames))

u.add_TopologyAttr(mda.core.topologyattrs.Tempfactors(np.zeros(u.atoms.n_atoms)))

print(" save {} begin frame and {} end frame with {} step".format(args.begin,args.end,args.step))
with mda.Writer(args.outpdb, multiframe=True, bonds=None, n_atoms=real_atoms.n_atoms) as PDB:
	for iframe in range(args.begin,args.end,args.step):
		ts = u.trajectory[iframe]
		# beta-factor
		temp = np.empty(0)
		for i_residue in range(n_residues):
			for i_atoms in range(len(u.residues.indices[i_residue])):
				temp = np.append(temp,np.abs(dih[iframe][i_residue])/10.)
		u.atoms.tempfactors = temp
		# position wrap each residues
		temp = np.empty(0)
		for i_residue in range(n_residues):
			list_atom = u.residues.indices[i_residue]
			res_a = mda.core.groups.AtomGroup(list_atom, u)
			cog = res_a.center_of_geometry()
			move_box = np.around(cog/u.dimensions[0:3]-0.5)
			#print(cog, move_box)
			temp = np.append(temp,u.atoms.positions[list_atom] - u.dimensions[0:3]*move_box)
			#print(temp)
		#print(temp.shape)
		u.atoms.positions = temp.reshape(-1,3)
		# for whole system
		# u.atoms.positions = u.atoms.positions - u.dimensions[0:3]*np.around(u.atoms.positions/u.dimensions[0:3]-0.5) 
		#print(real_atoms.tempfactors)
		#print(u.atoms.tempfactors)
		PDB.write(real_atoms.atoms)
		print(" ... current {} frame ...".format(ts.frame))
		
