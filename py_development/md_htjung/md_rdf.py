#!/usr/bin/env python2
# ver 0.1 - coding python by Hyuntae Jung on 7/19/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculate rdf only/except intramolecular pairs')
## args
parser.add_argument('-i', '--input', default='md.dcd', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='md.pdb', nargs='?', 
	help='pdb file for structure')
parser.add_argument('-step', '--step', default=1, nargs='?', type=int,
	help='iframe interval between frames you want to save')
parser.add_argument('-rmax', '--rmax', default=2, nargs='?', type=float,
	help='max distance of pairs to count, nm')
parser.add_argument('-nbin', '--nbin', default=100, nargs='?', type=int,
	help='rdf bin number')
parser.add_argument('-inc', '--include', default='except', nargs='?', 
	help='only including/except intramolecular pairs (except/only)')
parser.add_argument('-s1', '--sel1', nargs='?',
	help='file1 to select the atoms')
parser.add_argument('-s2', '--sel2', nargs='?',
	help='file2 to select the atoms')
parser.add_argument('-o', '--output', default='rdf_types', nargs='?', 
	help='output prefix')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

## import modules
import mdtraj as md
import numpy as np
import math

# load trajectory
traj = md.load(args.input,top=args.structure,stride=args.step)
top = traj.topology 
n_frames = len(traj.xyz)
print(" total frames = {}".format(len(traj.xyz)))

# read selection files
#  C - X : C for sel1, X for sel2
n_sel_command = []
for filename in [args.sel1, args.sel2]:
	try:
		open_file = open(filename, 'r')
	except IOError:
		raise IOError(" problem with opening ",filename)
	text1 = open_file.readline().strip()
	print(" select written in {}: \n  {} for selection command".format(filename,text1))
	n_sel_command.append(str(text1))
	open_file.close()

# make pair list
list_pairs = []
sel1 = top.select(n_sel_command[0])
sel2 = top.select(n_sel_command[1])
n_atom1 = len(sel1)
n_atom2 = len(sel2)
if args.include == 'except':
	for i_atom in sel1:
		for j_atom in sel2:
			if top.atom(i_atom).residue == top.atom(j_atom).residue:
				continue
			else:
				list_pairs.append((i_atom,j_atom))
elif args.include == 'only':
	for i_atom in sel1:
		for j_atom in sel2:
			if top.atom(i_atom).residue == top.atom(j_atom).residue:
				list_pairs.append((i_atom,j_atom))
			else:
				continue
else:
	raise ValueError(" wrong argument in args.include {}".format(str(args.include)))

list_pairs = np.array(list_pairs)
n_pairs = len(list_pairs)
print(" total # pairs = {} ".format(n_pairs))

# calculate distances 
dist = (md.compute_distances(traj,list_pairs)).flatten()

#print("remove elements where distance is greater than {}".format(args.rmax))
#del_list = np.where(dist >= args.rmax)
#print(del_list)
#print("{} pairs were removed.".format(len(del_list[0])))
#dist_red = np.delete(dist,del_list)
#print("max value {}, min value {} in the distance list".format(np.amax(dist_red), np.amin(dist_red)))
dist_red = dist
#print(len(dist_red))

# histogram
counts, edge_r = np.histogram(dist_red,bins=args.nbin,range=[0.0,args.rmax])
#print(counts)
#print(np.sum(counts))

# volume in each radial shell
vol = np.power(edge_r[1:],3) - np.power(edge_r[:-1],3)
vol *= 4/3.0 * np.pi
# Average number density
box_vol = np.average(traj.unitcell_volumes)
#print(box_vol)
density = len(dist_red) / box_vol

rdf = ( counts / density ) / vol
#print(rdf)

dist_x = 0.5*(edge_r[1:] + edge_r[:-1])
data = np.column_stack((dist_x,rdf))

np.savetxt(args.output, data, 
	header='x = distance [{},{}]' \
	.format(0,args.rmax), fmt='%f', comments='# ')
np.save(args.output, data)


# reference with MDTraj
#ref_data = md.compute_rdf(traj,pairs=list_pairs,r_range=(0.0,args.rmax),n_bins=args.nbin)
#print(ref_data[1])
#np.savetxt(args.output+'.ref', np.transpose(ref_data), 
#	header='x = distance [{},{}]' \
#	.format(0,args.rmax), fmt='%f', comments='# ')
