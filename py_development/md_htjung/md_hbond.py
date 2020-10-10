#!/usr/bin/env python2
# ver 0.1 - coding python by Hyuntae Jung on 05/24/2018
# ver 0.2 - updated for anion dihedral angles by Hyuntae Jung on 7/21/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculate 2D hbonds of ionic liquid in terms of dihedral angles of anions')
## args
parser.add_argument('-i', '--input', default='md.dcd', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='md.pdb', nargs='?', 
	help='pdb file for structure')
parser.add_argument('-stride', '--stride', default=100, nargs='?', type=int,
	help='only read n-th frames')
parser.add_argument('-dih', '--dih', default='No', nargs='?', 
	help='partial histogram in terms of anion dihedral angles (YES/NO)')
parser.add_argument('-sel', '--select', default='def_hbond.txt', nargs='?',
	help='command file to select an atom1 (of cation), atom2 (hydrogen bonded with atom1), and atom3 (electron donor in anion)')
parser.add_argument('-nbin_r', '--nbr', default=100, nargs='?', type=int,
	help='nbin of distances')
parser.add_argument('-rmax', '--rm', default=2.0, nargs='?', type=float,
	help='max of distances to draw')
parser.add_argument('-nbin_a', '--nba', default=360, nargs='?', type=int,
	help='nbin of angle')
parser.add_argument('-o', '--output', default='traj', nargs='?', 
	help='output prefix')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

args.output = args.output + '.rdf_2d'

## import modules
import mdtraj as md
import numpy as np
import math

# read trajectory
traj = md.load(args.input,top=args.structure,stride=args.stride)
top = traj.topology 
n_frames = len(traj.xyz)
print(" total frames = {}".format(len(traj.xyz)))

# catergorize dihedral angles of anions
if 'YES' in args.dih:
	# calculate dihedral angles of anions
	s1 = top.select("name Ctf")
	s2 = top.select("name Stf")
	s3 = top.select("name Stf1")
	s4 = top.select("name Ctf1")
	select_re = np.column_stack((s1,s2,s3,s4))
	dih_ani = md.compute_dihedrals(traj,select_re)
	# catergorize
	#  cisoid = -1, transoid = 0, intermediate = 1
	dih_ani = np.abs(dih_ani)*180.0/math.pi 
	dih_ani = np.where(dih_ani < 60, -1, dih_ani) 
	dih_ani = np.where(dih_ani > 120, 0, dih_ani)
	dih_ani = np.where(dih_ani > 0, 1, dih_ani)
	dih_ani = np.int_(dih_ani)
	n_ani = len(dih_ani[0]) # number of anions 
	
# read selection files
#  C - H ---- X : C for sel1, H for sel2, and X for sel3
n_sel_command = []
try:
	open_file = open(args.select, 'r')
except IOError:
	raise IOError(" problem with opening ",args.select)
text1 = open_file.readline().strip()
n_sel_command.append(str(text1))
text2 = open_file.readline().strip()
n_sel_command.append(str(text2))
text3 = open_file.readline().strip()
n_sel_command.append(str(text3))
print(" select written in {} for selection command: \n  {} \n {} \n {}".format(args.select,text1,text2,text3))
n_sel_command.append(str(text2))
open_file.close()

# make angle list
list_angles = []
list_angles_tran = []
list_angles_cis = []
list_angles_int = []
sel1 = top.select(n_sel_command[0])
sel2 = top.select(n_sel_command[1])
sel3 = top.select(n_sel_command[2])
n_atom1 = len(sel1)
n_atom2 = len(sel2)
n_atom3 = len(sel3)
if n_atom2%n_atom1 != 0:
	raise ValueError(" not consistent with # of atoms1 and atoms2")

n_cation = n_ani # assume n_cation = n_anion
n_dup1 = n_atom1/n_cation # #atom1s in single cation 
n_hs=int(n_atom2/n_atom1) # #atom2s bonded to single atom1
n_dup3 = n_atom3/n_ani # #atom3s in single anion

for i_atom in range(n_atom1):
	for k_atom in sel3:
		for repeat in range(n_hs):
			list_angles.append((sel1[i_atom],sel2[i_atom*n_hs+repeat],k_atom))

list_angles = np.array(list_angles)
n_angles = len(list_angles)
print(" total # angles = {} ".format(n_angles))

## make dihedral angle list which is the same size as angl
if 'YES' in args.dih:
	tmp_dih = np.repeat(dih_ani,n_hs,axis=1)
	tmp_dih = np.repeat(tmp_dih,n_dup3,axis=1)
	tmp_dih = np.repeat(tmp_dih,n_atom1,axis=0)
	dih_2_angl = tmp_dih.flatten()
	# get indices for trans, cis, and intermediates
	ind_cis = np.where(dih_2_angl < -0.5)
	ind_tran = np.where(dih_2_angl > 0.5)
	ind_intm = np.where(dih_2_angl == 0)

## calculate angles and distances 
angl = (md.compute_angles(traj,list_angles)).flatten()
if len(angl) != len(dih_2_angl):
	raise ValueError("different size of angl and dih_2_angl")

# recalculate NaN values of angles
import copy
for nan_item in np.argwhere(np.isnan(angl)).reshape(-1):
	i_frame = int(nan_item/n_angles)
	i_angle = nan_item%n_angles
	print(" Nan at {} th frame and {} th angle".format(i_frame,i_angle))
	print(" <-- {} th atoms".format(list_angles[i_angle]))
	i_abc = list_angles[i_angle]
	a = traj.xyz[i_frame][i_abc[0]]
	b = traj.xyz[i_frame][i_abc[1]]
	c = traj.xyz[i_frame][i_abc[2]]
	print(" <-- position: \n {} \n {} \n {}".format(a,b,c))
	ba = a - b
	bc = c - b
	cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
	print(" distance= {}".format(np.linalg.norm(bc)))
	angle = np.arccos(cosine_angle)
	print(" get correct value from NaN, {} (rad) {} (deg)".format(angle, angle*180.0/np.pi))
	angl[nan_item] = copy.copy(angle)

# calculate distance
dist = (md.compute_distances(traj,list_angles[:,[1,2]])).flatten()

## histogram
counts_2d, edge_r, edge_a = np.histogram2d(dist,angl,bins=[args.nbr,args.nba],range=[[0.0,args.rm],[0.0,np.pi]])
# volume in each radial shell
vol = np.power(edge_r[1:],3) - np.power(edge_r[:-1],3)
vol *= 4/3.0 * np.pi
# Average number density
box_vol = np.average(traj.unitcell_volumes)
density = len(dist) / box_vol
rdf_2d = ( counts_2d*args.nba / n_frames ) / vol[:,None]
np.savetxt(args.output, np.transpose(rdf_2d), 
	header='x = distance [{},{}], y = abs angle [{},{}]' \
	.format(0,args.rm,0,np.pi), fmt='%f', comments='# ')
np.save(args.output, np.transpose(rdf_2d))

if 'YES' in args.dih:
	for list_indice, prefix_file in zip([ind_cis, ind_intm, ind_tran], [".cis", ".intm", ".tran"]):
		counts_2d, edge_r, edge_a = np.histogram2d(dist[list_indice],angl[list_indice],bins=[args.nbr,args.nba],range=[[0.0,args.rm],[0.0,np.pi]])

		# volume in each radial shell
		vol = np.power(edge_r[1:],3) - np.power(edge_r[:-1],3)
		vol *= 4/3.0 * np.pi
		# Average number density
		density = len(dist[list_indice]) / box_vol
		rdf_2d = ( counts_2d*args.nba / density ) / vol[:,None]
		np.savetxt(args.output+prefix_file, np.transpose(rdf_2d), 
			header='x = distance [{},{}], y = abs angle [{},{}]' \
			.format(0,args.rm,0,np.pi), fmt='%f', comments='# ')
		np.save(args.output+prefix_file, np.transpose(rdf_2d))

