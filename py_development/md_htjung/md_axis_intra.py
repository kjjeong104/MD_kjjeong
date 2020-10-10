#!/usr/bin/env python2
# ver 0.1 - coding python by Hyuntae Jung on 7/22/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculate RADFs of molecular vectors (or axis)')
## args
parser.add_argument('-i', '--input', default='md.dcd', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='md.pdb', nargs='?', 
	help='pdb file for structure')
parser.add_argument('-step', '--step', default=1, nargs='?', type=int,
	help='iframe interval between frames you want to save')
parser.add_argument('-rmax', '--rmax', default=2, nargs='?', type=float,
	help='max distance of pairs to count, nm')
parser.add_argument('-nbin_r', '--nbr', default=100, nargs='?', type=int,
	help='rdf bin number')
parser.add_argument('-nbin_a', '--nba', default=360, nargs='?', type=int,
	help='nbin of angle')
parser.add_argument('-sel', '--select', nargs='?',
	help='file to select the two atoms of molecule1 and molecule2')
parser.add_argument('-dih', '--dih', default="FF", nargs='?',
	help='partial RADFs with dihedral angles of anion (TT/FT/FF); pRADFs applied to sel1 and sel1? T/F = True/False')
parser.add_argument('-o', '--output', default='axis_rdf', nargs='?', 
	help='output prefix')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

args.output = args.output + '.axis_rdf_2d'

## import modules
import mdtraj as md
import numpy as np
import math
# track RuntimeWarning 
import warnings 
warnings.simplefilter('error')

# read trajectory
traj = md.load(args.input,top=args.structure,stride=args.step)
top = traj.topology 
n_frames = len(traj.xyz)
print(" total # loaded frames = {}".format(len(traj.xyz)))

# catergorize and calculate dihedral angles of anions
if 'T' in args.dih:
	s1 = top.select("name Ctf")
	s2 = top.select("name Stf")
	s3 = top.select("name Stf1")
	s4 = top.select("name Ctf1")
	select_re = np.column_stack((s1,s2,s3,s4))
	dih_ani = md.compute_dihedrals(traj,select_re)
	# catergorize
	#  cisoid = -1, transoid = 0, intermediate = 2
	dih_ani = np.abs(dih_ani)*180.0/math.pi 
	dih_ani = np.where(dih_ani < 60, -1, dih_ani) 
	dih_ani = np.where(dih_ani > 120, 0, dih_ani)
	dih_ani = np.where(dih_ani > 0, 2, dih_ani)
	dih_ani = np.int_(dih_ani)
	n_ani = len(dih_ani[0]) # number of anions 
	if args.dih == 'TT':
		print("pRADFs for both sel1 and sel2")
	elif args.dih == 'FT':
		print("pRADFs only for sel2")
	else:
		raise ValueError(" wrong arguments in args.dih, {}".format(args.dih))
elif args.dih == 'FF':
	print("No pRADFs")
else:
	raise ValueError(" wrong arguments on args.dih, {}".format(args.dih))

# read selection files
#  A-B and C-D for sel1, and sel2
n_sel_command = []
try:
	open_file = open(args.select, 'r')
except IOError:
	raise IOError(" problem with opening ",args.select)
text1 = open_file.readline().strip()
n_sel_command.append(str(text1))
text2 = open_file.readline().strip()
n_sel_command.append(str(text2))
print(" select written in {} for selection command: \n  {} \n  {}".format(args.select,text1,text2))
open_file.close()

# make pair list
list_pairs_sel1 = []
list_pairs_sel2 = []
sel1 = top.select(n_sel_command[0])
sel2 = top.select(n_sel_command[1])
n_atom1 = len(sel1)
n_atom2 = len(sel2)
if n_atom1%2 != 0:
	raise ValueError(" wrong selection (not pairs) for select 1")

if n_atom2%2 != 0:
	raise ValueError(" wrong selection (not pairs) for select 2")

list_vector_pairs = []
if n_sel_command[0] == n_sel_command[1]:
	for i in range(int(n_atom1/2)):
		list_pairs_sel1.append((sel1[2*i],sel1[2*i+1]))
	n_pairs = len(list_pairs_sel1)
	for i_pair in range(n_pairs):
		for j_pair in range(i_pair+1,n_pairs):
			list_vector_pairs.append((list_pairs_sel1[i_pair],list_pairs_sel1[j_pair]))

	## get indices for each combinations
	if 'TT' in args.dih: # anion-anion case; trans, cis, interm combinations
		dih_comb=[]
		for iframe in dih_ani:
			n_mol = len(iframe)
			for i in range(n_mol):
				for ii in (iframe[i]+iframe[i+1:n_mol]):
					dih_comb.append(ii)

		dih_comb = np.array(dih_comb)
		print(dih_comb,dih_comb.shape)
		t_t = np.where((dih_comb > -0.5) & (dih_comb < 0.5))  # trans-trans = 0
		t_c = np.where((dih_comb > -1.5) & (dih_comb < -0.5)) # trans-cis = -1
		t_i = np.where((dih_comb > 1.5) & (dih_comb < 2.5))  # trans-interm = 2
		c_c = np.where(dih_comb < -1.5) # cis-cis = -2
		print(c_c)
		c_i = np.where((dih_comb > 0.5) & (dih_comb < 1.5))  # cis-interm = 1
		i_i = np.where(dih_comb > 3.5)  # interm-interm = 4
		print("#trans-trans, #trans-cis, #trans-interm, #cis-cis, #cis-interm, #interm-interm \n {} {} {} {} {} {}"
			.format(len(t_t[0]),len(t_c[0]),len(t_i[0]),len(c_c[0]),len(c_i[0]),len(i_i[0])))		
else:
	for i in range(int(n_atom1/2)):
		list_pairs_sel1.append((sel1[2*i],sel1[2*i+1]))

	for i in range(int(n_atom2/2)):
		list_pairs_sel2.append((sel2[2*i],sel2[2*i+1]))

	for i_pair in list_pairs_sel1:
		for j_pair in list_pairs_sel2:
			if i_pair == j_pair:
				raise ValueError(" somthing wrong with some identical pairs found")
			else:
				list_vector_pairs.append((i_pair,j_pair))

	if 'FT' in args.dih: # cation-anion case; trans, cis, interm for anion
		n_list_pairs_sel1 = len(list_pairs_sel1)
		dih_comb = (np.repeat(dih_ani,n_list_pairs_sel1,axis=0)).flatten()
		ca_t = np.where((dih_comb > -0.5) & (dih_comb < 0.5))  # any-trans = 0
		ca_c = np.where(dih_comb < -0.5) # any-cis = -1
		ca_i = np.where(dih_comb > 1.5)  # any-interm = 2
		print("#any-trans, #any-cis, #any-interm \n {} {} {}"
			.format(len(ca_t[0]),len(ca_c[0]),len(ca_i[0])))		

list_vector_pairs = np.array(list_vector_pairs)
n_vector_pairs = len(list_vector_pairs)
print(" total # vector pairs = {} ".format(n_vector_pairs))

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

def run_histrogram_2d(list_indice, prefix_file):
	## histogram
	if list_indice is not None:
		#print((dist[list_indice]).shape)
		counts_2d, edge_r, edge_a = np.histogram2d(dist[list_indice],angl[list_indice],
			bins=[args.nbr,args.nba],range=[[0.0,args.rmax],[0.0,np.pi]])
		#print(np.sum(counts_2d))
	else:
		counts_2d, edge_r, edge_a = np.histogram2d(dist,angl,
			bins=[args.nbr,args.nba],range=[[0.0,args.rmax],[0.0,np.pi]])
	
	# volume in each radial shell
	vol = np.power(edge_r[1:],3) - np.power(edge_r[:-1],3)
	vol *= 4/3.0 * np.pi
	
	# Average number density
	box_vol = np.average(traj.unitcell_volumes)
	if list_indice is not None:
		density = len(dist[list_indice]) / box_vol
	else:
		density = len(dist) / box_vol
	rdf_2d = ( counts_2d*args.nba / density ) / vol[:,None]
	if list_indice is not None:
		outfile=args.output+prefix_file
	else:
		outfile=args.output
	np.savetxt(outfile, np.transpose(rdf_2d), 
		header='x = distance [{},{}], y = abs angle [{},{}]' \
		.format(0,args.rmax,0,np.pi), fmt='%f', comments='# ')
	np.save(outfile, np.transpose(rdf_2d))

if args.dih == 'TT': # anion-anion case;
	for list_indice, prefix_file in zip([t_t, t_c, t_i, c_c, c_i, i_i], [".t_t", ".t_c", ".t_i", ".c_c", ".c_i", ".i_i"]):
		run_histrogram_2d(list_indice, prefix_file)
elif args.dih == 'FT': # cations (or total anions) - anion case
	for list_indice, prefix_file in zip([ca_t, ca_c, ca_i], [".ca_t", ".ca_c", ".ca_i"]):
		run_histrogram_2d(list_indice, prefix_file)
else:
	run_histrogram_2d(None, None)
