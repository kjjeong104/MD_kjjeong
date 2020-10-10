#!/usr/bin/env python2
# ver 0.1 - coding python by Hyuntae Jung on 12/29/2017

import argparse
parser = argparse.ArgumentParser(
	formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
	description='calculate dihedral angles')
## args
parser.add_argument('-i', '--input', default='md.dcd', nargs='?', 
	help='input trajectory file')
parser.add_argument('-s', '--structure', default='md.pdb', nargs='?', 
	help='pdb file for structure')
parser.add_argument('-b', '--begin', default=0, nargs='?', type=int,
	help='starting iframe')
parser.add_argument('-mol', '--mol', default='cation', nargs='?',
	help='select the molecule type (cation/anion)')
parser.add_argument('-o', '--output', default='traj', nargs='?', 
	help='output prefix')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

args.odih = args.output + '.dih'
args.odih_raw = args.output + '.raw.dih'

## import modules
import mdtraj as md
import numpy as np
import math

traj = md.load(args.input,top=args.structure)
top = traj.topology 
print(" total frames = {}".format(len(traj.xyz)))
if 'cation' in args.mol:
	#select = top.select("name N or name C3 or name C4 or name O or name C5")
	s1 = top.select("name N")
	s2 = top.select("name C3")
	s3 = top.select("name C4")
	s4 = top.select("name O or name C5")
	select_re = np.column_stack((s1,s2,s3,s4))
elif 'anion' in args.mol:
	#select = top.select("name Ctf or name Stf or name Stf1 or name Ctf1")
	s1 = top.select("name Ctf")
	s2 = top.select("name Stf")
	s3 = top.select("name Stf1")
	s4 = top.select("name Ctf1")
	select_re = np.column_stack((s1,s2,s3,s4))
else:
	raise ValueError("arugment {} is not correct (cation/anion)".format(args.mol))
print(" {} torions you select".format(len(select_re)))
dih = md.compute_dihedrals(traj,select_re)
# cut frames
print(" starting with {} th frames".format(args.begin))
dih_resize = dih[args.begin:]
# histogram
prop_dih, radian_dih = np.histogram(np.abs(dih_resize),range=(0,math.pi),bins=100)
n_dih = np.sum(prop_dih)
prop_dih = prop_dih/float(n_dih)
deg_dih = radian_dih[0:-1]*180.0/math.pi
data = np.column_stack((deg_dih,prop_dih))

np.savetxt(args.odih, data, 
	header='dihedral angles of {} from {} th frames' \
	.format(args.mol,args.begin), fmt='%f', comments='# ')
np.save(args.odih, data)
np.save(args.odih_raw, dih)
