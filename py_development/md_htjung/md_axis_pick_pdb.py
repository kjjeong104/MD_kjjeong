#!/usr/bin/env python2
# ver 0.1 - coding python from md_axis.py by Hyuntae Jung on 8/6/2018

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='pick up information from RADFs of molecular vectors')
## args
parser.add_argument('-i', '--input', default='md.dcd', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='md.pdb', nargs='?', 
	help='pdb file for structure')
parser.add_argument('-step', '--step', default=1, nargs='?', type=int,
	help='iframe interval between frames you want to save')
parser.add_argument('-sel', '--select', nargs='?',
	help='file to select the two atoms of molecule1 and molecule2')
parser.add_argument('-dih', '--dih', default="FF", nargs='?',
	help='partial RADFs with dihedral angles of anion (TT/FT/FF); pRADFs applied to sel1 and sel1? T/F = True/False')
parser.add_argument('-pd', '--pd', default=0.43, nargs='?', type=float,
	help='average distance of pairs to pick, nm [value-pdi:value+pdi]')
parser.add_argument('-pdi', '--pdi', default=0.03, nargs='?', type=float,
	help='distance interval range')
parser.add_argument('-pa', '--pa', default=152, nargs='?', type=float,
	help='average distance of pairs to pick, nm [value-10:value+10]')
parser.add_argument('-pai', '--pai', default=5, nargs='?', type=float,
	help='angle interval range')
parser.add_argument('-o', '--output', default='axis_rdf', nargs='?', 
	help='output prefix')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

args.output = args.output + '.axis_rdf_2d.pick'

## import modules
import mdtraj as md
import numpy as np
import math
import re
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
dist_2d = compute_vector_distances(traj,list_vector_pairs)
print(" done with computing vector distances")
angl_2d = compute_vector_angles(traj,list_vector_pairs)
print(" done with computing vector angles")

# find what we want to pick up
dist_2d_pick = (dist_2d > args.pd-args.pdi) & (dist_2d < args.pd+args.pdi)
lst_angl = (args.pa-args.pai)*math.pi/180.0
grt_angl = (args.pa+args.pai)*math.pi/180.0
angl_2d_pick = (angl_2d > lst_angl) & (angl_2d < grt_angl)
data_2d_pick = np.where(dist_2d_pick & angl_2d_pick)

# info output
n_items = len(data_2d_pick[0])
for i in range(n_items):
	print("we found:")
	i_pair = data_2d_pick[0][i]
	print(" {} th pair".format(i_pair))
	(a,b),(c,d) = list_vector_pairs[i_pair]
	text1 = re.split(r'[/s-]+', str(top.atom(a)))
	text2 = re.split(r'[/s-]+', str(top.atom(c)))
	resname1 = text1[0][0:3]
	resid1 = text1[0][3:]
	resname2 = text2[0][0:3]
	resid2 = text2[0][3:]
	print(" residues = {}th {} - {}th {} in vmd".format(resid1,resname1,resid2,resname2))
	if resname2 == 'Tf2':
		resid1 = int(resid1)
		resid2 = int(resid2) + 200
	print(" residues = {}th {} - {}th {} in mdtraj".format(resid1-1,resname1,resid2-1,resname2))
	sel_t1 = "(resname "+str(resname1)+" and resid "+str(resid1-1)+")"
	sel_t2 = "(resname "+str(resname2)+" and resid "+str(resid2-1)+")"
	pick_sel = top.select("("+str(sel_t1)+" or "+str(sel_t2)+") and not (name =~ 'D.*')")
	traj_slice = traj.atom_slice(pick_sel)
	
	i_frame = data_2d_pick[1][i]
	print(" {} th frame (after stride {} frames)".format(i_frame,args.step))
	print(" xyz for {} th atom of {} residue: \n {}".format(a,top.atom(a),traj.xyz[i_frame,a]))
	out_traj = traj_slice[i_frame]
	
	# closed form by pbc
	origin = out_traj.xyz[0][0]
	out_traj.xyz[0] = out_traj.xyz[0] - origin
	n_atoms_sel = len(out_traj.xyz[0])
	for i_pos in range(n_atoms_sel) :
		out_traj.xyz[0][i_pos][:] = out_traj.xyz[0][i_pos][:] - np.around(out_traj.xyz[0][i_pos][:]/out_traj.unitcell_lengths[0][:])*out_traj.unitcell_lengths[0][:]
	out_traj.save_gro("pick.gro")
	traj[i_frame].save_gro("pick_tot.gro")
	print(" saved in pick.gro and pick_tot.gro")
	input("Press Enter to continue...")
	

print("Done")
