#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 05/11/2018
################
## Transport and Helfand moments in the Lennard-Jones fluid. I. Shear viscosity
## S. Viscardy, J. Servantie, and P. Gaspard
##  Citation: The Journal of Chemical Physics 126, 184512 (2007); doi: 10.1063/1.2724820
##  View online: https://doi.org/10.1063/1.2724820
################

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculation shear viscosity cubic system using Einstein relation by Helfand')
## args
parser.add_argument('-i', '--input', default='traj.trr', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='topol.tpr', nargs='?', 
	help='.tpr or .gro structure file')
parser.add_argument('-temp', '--temp', default=298, nargs='?', type=float,
	help='temperature (unit = kelvin)')
parser.add_argument('-step', '--stepsize', default=1, nargs='?', type=float,
	help='time interval between neighbor frames in trajectory (unit = ps)')
parser.add_argument('-o', '--output', default='vis', nargs='?', 
	help='output prefix for instantaneous viscosity')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
# read args
args = parser.parse_args()
# check args
print(" input arguments: {0}".format(args))

## import modules
import sys
sys.path.append('/home/htjung/Utility/python')
import hjung
from hjung import *
import numpy as np
import MDAnalysis as md
import copy

args.output = args.output + '.out'

## timer
start_proc, start_prof = hjung.time.init()

# read trajectory then calculate center of mass
#  assume that the positions are not jumped by periodic boundary condition
u = md.Universe(args.structure, args.input)
n_frames = len(u.trajectory)
n_atoms = len(u.atoms)
vol = u.dimensions[0]*u.dimensions[1]*u.dimensions[2] # assume NVT ensemble

# generate edited mass array 
mass_array = copy.copy(u.atoms.masses)
for i in range(n_atoms):
	if u.atoms.types[i] is '': # drude particle 
		mass_array[i] = 0.40
	elif u.atoms.types[i] is 'D': # drude particle
		mass_array[i] = 0.40
	elif u.atoms.types[i] is 'H': # hydrogen
		continue
	else:
		mass_array[i] = mass_array[i] - 0.40

# calculate pos and vel for center of mass every residue
n_residues = len(u.residues)
pos = np.empty((n_frames,n_residues,3))
vel = np.empty((n_frames,n_residues,3))
tmass = np.empty(n_residues)
mod_frame = hjung.time.process_init()
for i_frame in range(n_frames):
	ts = u.trajectory[i_frame]
	for i_residue in range(n_residues):
		#print(" ... residue {} ... ".format(i_residue))
		select_residue = u.select_atoms("resid "+str(i_residue+1)+" and (not (name =~ 'D'))") # select a residue exclude drude particles
		n_atoms_residue = len(select_residue)
		if n_atoms_residue == 0:
			raise ValueError("wrong argument resid because of no selection. Check topology file (resid)")
		#else:
		#	print(" {} atoms for residue {}".format(n_atoms_residue,i_residue))
		# calculate center of mass
		if i_frame == 0:
			total_mass = np.sum(mass_array[select_residue.indices])
			tmass[i_residue] = total_mass
		pos[i_frame,i_residue] = np.einsum('i,ij->j',mass_array[select_residue.indices],select_residue.positions)/tmass[i_residue]
		vel[i_frame,i_residue] = np.einsum('i,ij->j',mass_array[select_residue.indices],select_residue.velocities)/tmass[i_residue]
	mod_frame = hjung.time.process_print(i_frame+1,n_frames,mod_frame)

# recenter to center of mass of residues
total_mass = np.sum(tmass)
for i_frame in range(n_frames):
	pos_c = np.einsum('i,ij->j',tmass,pos[i_frame])/total_mass
	vel_c = np.einsum('i,ij->j',tmass,vel[i_frame])/total_mass
	#print(pos_c,vel_c)
	pos[i_frame] = pos[i_frame] - pos_c
	vel[i_frame] = vel[i_frame] - vel_c

# generate dt list for viscosity
list_dt = []
basic_list = np.arange(1,10)
while True: 
	if basic_list[0] >= n_frames:
		break
	for element in basic_list:
		if element < n_frames:
			list_dt.append(element)
		else:
			break
	basic_list = basic_list*10
list_dt = np.array(list_dt[:-1])
n_dt = len(list_dt)

# calculate large parathesis in eq 17. in paper above
moments = pos[:,:,0]*vel[:,:,1]*tmass # 2d

#from timeit import default_timer as timer
## method1: calc. all time delays, but fast algorithm
##  speed= 0.037 sec per 200 frames
#start = timer()
diff_moments, idx = hjung.analyze.einstein_time_relation(moments,mode='idx')
momentum = np.empty((n_dt,2))
j=0
for i in list_dt:
	temp = np.sum(diff_moments[idx[:-1-i]+i-1],axis=1)**2
	#print(temp)
	#print(vol)
	momentum[j,0] = np.average(temp)/float(i)
	momentum[j,1] = np.std(temp)/float(i)
	j = j + 1
#end = timer()
#print(end - start)

## method2: calc only list_dt 
##  speed= 0.042 sec per 200 frames
#start = timer()
#j = 0
#momentum = np.empty((n_dt,2))
#for idt in list_dt:
#	n_origin = len(moments)-idt
#	temp = np.empty(n_origin)
#	for origin_t in range(n_origin):
#		temp[origin_t] = np.sum(moments[idt+origin_t] - moments[origin_t])**2
#	momentum[j,0] = np.average(temp)/float(idt)
#	momentum[j,1] = np.std(temp)/float(idt)
#	j = j + 1
#end = timer()
#print(end - start)

momentum = momentum/(2.*vol*6.0221417930*8.3144598*args.temp*args.stepsize) # unit = centipoise
# later it should divide by {} picoseconds (for one frame time size in trajecotry file)

data = np.column_stack((list_dt,momentum))
print(data)
np.savetxt(args.output,data)
np.save(args.output,data)

## timer
hjung.time.end_print(start_proc, start_prof)