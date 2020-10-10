#!/usr/bin/env python2
# ver 0.1 - coding python from md_axis.py by Hyuntae Jung on 8/15/2018 

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculate (self) rotational correlation of molecular vectors (or axis)')
## args
parser.add_argument('-i', '--input', default='md.dcd', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='md.pdb', nargs='?', 
	help='pdb file for structure')
parser.add_argument('-step', '--step', default=1, nargs='?', type=int,
	help='iframe interval between frames you want to save')
parser.add_argument('-ion', '--ion', default='chol', nargs='?',
	help='what ion to analyze (chol/tfsi)')
parser.add_argument('-sel', '--select', nargs='?',
	help='file to select two atoms for each ion')
parser.add_argument('-o', '--output', default='axis_rot', nargs='?', 
	help='output prefix')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

args.output = args.output + '.axis_rot'

## import modules
import mdtraj as md
import numpy as np
import math
# track RuntimeWarning 
import warnings 
warnings.simplefilter('error')
# load hjung module
sys.path.append('/home/htjung/Utility/python/')
import hjung
from hjung import *

# read trajectory
traj = md.load(args.input,top=args.structure,stride=args.step)
top = traj.topology 
n_frames = len(traj.xyz)
print(" total # loaded frames = {}".format(len(traj.xyz)))

# read selection files
#  A-B atom list (for one line)
n_sel_command = hjung.io.read_text_line(args.select,2)

# make pair list
sel1 = top.select(n_sel_command[0])
n_atom1 = len(sel1)

list_vector_pairs = []
for i in range(int(n_atom1/2)):
	list_vector_pairs.append((sel1[2*i],sel1[2*i+1]))
list_vector_pairs = np.array(list_vector_pairs)
n_vector_pairs = len(list_vector_pairs)
print(" total # vector pairs = {} ".format(n_vector_pairs))


## get com of resid
for i,j in list_vector_pairs:
	#show = process_init()
	print(" ... resid {} ... ".format(resid_index))
	select_atoms = top.select(resid+str(resid_index)+" and (not (name =~ 'D'))") # select a residue exclude drude particles
	# check number of selected atoms
	if len(select_atoms) == 0:
		raise ValueError("wrong argument resid because of no selection")
	if natoms != len(select_atoms):
		if natoms != 0:
			raise RuntimeError("number atoms of resid {} is changed during selection. Probably wrong arugment resid range".format(resid_index))
		else:
			natoms = len(select_atoms)
	traj_resid_wo_drude = md.load(args.input,top=args.structure,atom_indices=select_atoms) # load trajectory
	com = md.compute_center_of_mass(traj_resid_wo_drude) # calc com
	com = com[1:]
	delta_com = com - np.mean(com, axis=0)
	for i in range(3):
		acf_com = ndimage.correlate(delta_com[:,i],delta_com[:,i],mode='constant')
		acf_data_1d /= (com[:,i].var()*len(delta_com[:,i])) # normalize
	#print(com[0])
	resid_msd, resid_msd_iso = calc_msd_com(com,basic_list) # calc msd
	# sum values
	msd = msd + resid_msd
	msd_iso = msd_iso + resid_msd_iso
	msd_iso_std = msd_iso_std + resid_msd_iso**2
	# print
	#show = process_print(resid_index,args.resid_start+n_resid,show)


## compute vector distance
def compute_vector_distances(traj, pairs):
	n_pairs = len(pairs)
	n_frames = len(traj)
	out = np.zeros((n_pairs,n_frames))
	for i in range(n_pairs):
		if i%(int(n_pairs/10)) == 0:
			print(" ... {} th pair ... ".format(i))
		(a,b),(c,d) = pairs[i]
		axis1 = (traj.xyz[:,b] + traj.xyz[:,a])/2.0 # not necessary to consider pbc condition because a and b atoms are within single molecule
		axis2 = (traj.xyz[:,d] + traj.xyz[:,c])/2.0
		dist = axis2-axis1
		for xyz in range(0,2):
			dist[:,xyz] = dist[:,xyz] - np.around(dist[:,xyz]/traj.unitcell_lengths[:,xyz])*traj.unitcell_lengths[:,xyz]
		out[i] = np.linalg.norm(dist,axis=1)
	return out

def compute_vector_angles(traj, pairs):
	n_pairs = len(pairs)
	n_frames = len(traj)
	out = np.zeros((n_pairs,n_frames))
	for i in range(n_pairs):
		if i%(int(n_pairs/10)) == 0:
			print(" ... {} th pair ... ".format(i))
		(a,b),(c,d) = pairs[i]
		axis1 = traj.xyz[:,b] - traj.xyz[:,a] # not necessary to consider pbc condition because a and b atoms are within single molecule
		axis2 = traj.xyz[:,d] - traj.xyz[:,c]
		dot_axis = np.sum(axis1*axis2,axis=1)
		length_dot = np.linalg.norm(axis1,axis=1)*np.linalg.norm(axis2,axis=1)
		try:
			out[i] = np.arccos(dot_axis/length_dot)
		except RuntimeWarning: # it happens when valid digit error
			#print(i)
			#print(np.amax(dot_axis),length_dot[np.argmax(dot_axis)])
			#print(np.amin(length_dot),dot_axis[np.argmin(length_dot)])
			#print(np.amax(dot_axis/length_dot),dot_axis[np.argmax(dot_axis/length_dot)],length_dot[np.argmax(dot_axis/length_dot)])
			#print(np.amin(dot_axis/length_dot),dot_axis[np.argmin(dot_axis/length_dot)],length_dot[np.argmin(dot_axis/length_dot)])
			if np.amin(dot_axis/length_dot) < -1.0:
				print("RuntimeWarning1 corrected")
				out[i,np.argmin(dot_axis/length_dot)] = np.arccos(-1.0)
			if np.amax(dot_axis/length_dot) > 1.0:
				print("RuntimeWarning2 corrected")
				out[i,np.argmax(dot_axis/length_dot)] = np.arccos(1.0)
	return out

## calculate vector distance and angles
dist = (compute_vector_distances(traj,list_vector_pairs)).flatten()
print(" done with computing vector distances")
angl = (compute_vector_angles(traj,list_vector_pairs)).flatten()
print(" done with computing vector angles")
