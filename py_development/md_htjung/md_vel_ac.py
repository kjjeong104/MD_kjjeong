#!/usr/bin/env python2
# ver 0.1 - coding python by Hyuntae Jung on 7/29/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculate com velocities from velocity trajectory file')
## args
parser.add_argument('-i', '--input', default='md_vel.dcd', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='md.pdb', nargs='?', 
	help='pdb file for structure')
parser.add_argument('-resid_start', '--resid_start', nargs='?', type=int, 
	help='index of starting resid to analyze (starting 0)')
parser.add_argument('-resid_end', '--resid_end', nargs='?', type=int, 
	help='index of ending resid to analyze')
parser.add_argument('-o', '--output', default='traj', nargs='?', 
	help='output prefix for center of mass trajectory')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

args.ovel = args.output + '.com_vel'
args.ovel_iso = args.output + '.com_vel_iso'

## import modules
import mdtraj as md
import numpy as np
from scipy import ndimage
sys.path.append('/home/htjung/Utility/python/')
import hjung
from hjung import *

## make atom_indices for each molecules
n_resid = args.resid_end - args.resid_start + 1
traj = md.load(args.input,top=args.structure)
n_frames = len(traj)
top = traj.topology 
resid="resid "

## check
if args.resid_end >= top.n_residues:
	raise ValueError("argument resid_end is beyond the range of resid in trajectory")
if args.resid_start < 0:
	raise ValueError("argument resid_start is beyond the range of resid in trajectory")

com_vel_iso = np.zeros((n_frames+1,n_resid))
com_vel_x = np.zeros((n_frames+1,n_resid))
com_vel_y = np.zeros((n_frames+1,n_resid))
com_vel_z = np.zeros((n_frames+1,n_resid))

msd = np.zeros((dtn,3)) # msd for x, y, z direction 
msd_iso = np.zeros(dtn) # msd for square distance (isotropy)
msd_iso_std = np.zeros(dtn) # msd_std for square distance (isotropy)
natoms = 0 # check any selection error later

## calculate eigen vectors

## get com of resid
for resid_index in range(args.resid_start,args.resid_start+n_resid):
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

## calc avg, std of msd
msd_avg = msd/float(n_resid)
msd_iso_avg = msd_iso/float(n_resid)
msd_iso_std = msd_iso_std/float(n_resid)
msd_iso_std = np.sqrt(msd_iso_std - msd_iso_avg**2)

#basic_list = np.array(basic_list)
data_msd_iso = np.stack((basic_list,msd_iso_avg,msd_iso_std),axis=1)
data_msd = np.column_stack((basic_list,msd_avg))

np.savetxt(args.omsd, data_msd, 
	header='msd (x,y,z) with selection resid {} to {}' \
	.format(args.resid_start,args.resid_end), fmt='%f', comments='# ')
np.save(args.omsd, data_msd)
np.savetxt(args.omsd_iso, data_msd_iso, 
	header='msd (iso) with selection resid {} to {}' \
	.format(args.resid_start,args.resid_end), fmt='%f', comments='# ')
np.save(args.omsd_iso, data_msd_iso)
