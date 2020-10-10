#!/usr/bin/env python3

# intermolecular contact (vdw fraction) calculator
# input : concatenated xyz file(dimers), number of atoms for monomer 1 and  monomer 2
# output : column file tells about minimum vdw fraction
# calculate vdw fraction of all intermolecular atomic contact, and give minimum value

# These van der waals radii are taken from
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
# which references 
# A. Bondi (1964). "van der Waals Volumes and Radii". J. Phys. Chem. 68: 441. doi:10.1021/j100785a001 and doi:10.1021/jp8111556. 
# M. Mantina et al. (2009). "Consistent van der Waals Radii for the Whole Main Group". J. Phys. Chem. A. 113 (19): 5806--12. doi:10.1021/jp8111556
# Where no van der Waals value is known, a default of 2 angstroms is used.
# However, because certain atoms in biophysical simulations have a high 
# chance of being completely ionized, we have decided to give the 
# following atoms their ionic radii their ionic radii:
# +2: Be, Mg, Ca, Ba
# +1: Li, Na, K, Cs
# -1: Cl
# These ionic radii are were taken from:
# Shannon, R. D. Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. Acta Crystallographica Section A 32, 751--767 (1976). doi:10.1107/S0567739476001551
# For most atoms, adding electrons usually doesn't change the radius much 
# (<10%), while removing them changes it substantially (>50%). Further, 
# when atoms like N, S, and P, are positive, they are bound to atoms in such 
# a way that would "hide" their radii anyway. We have therefore chosen to just 
# use their vdW radii.

_ATOMIC_RADII = {'H'   : 0.120, 'He'  : 0.140, 'Li'  : 0.076, 'Be' : 0.059, 
                 'B'   : 0.192, 'C'   : 0.170, 'N'   : 0.155, 'O'  : 0.152, 
                 'F'   : 0.147, 'Ne'  : 0.154, 'Na'  : 0.102, 'Mg' : 0.086, 
                 'Al'  : 0.184, 'Si'  : 0.210, 'P'   : 0.180, 'S'  : 0.180, 
                 'Cl'  : 0.181, 'Ar'  : 0.188, 'K'   : 0.138, 'Ca' : 0.114, 
                 'Sc'  : 0.211, 'Ti'  : 0.200, 'V'   : 0.200, 'Cr' : 0.200, 
                 'Mn'  : 0.200, 'Fe'  : 0.200, 'Co'  : 0.200, 'Ni' : 0.163, 
                 'Cu'  : 0.140, 'Zn'  : 0.139, 'Ga'  : 0.187, 'Ge' : 0.211, 
                 'As'  : 0.185, 'Se'  : 0.190, 'Br'  : 0.185, 'Kr' : 0.202, 
                 'Rb'  : 0.303, 'Sr'  : 0.249, 'Y'   : 0.200, 'Zr' : 0.200, 
                 'Nb'  : 0.200, 'Mo'  : 0.200, 'Tc'  : 0.200, 'Ru' : 0.200, 
                 'Rh'  : 0.200, 'Pd'  : 0.163, 'Ag'  : 0.172, 'Cd' : 0.158, 
                 'In'  : 0.193, 'Sn'  : 0.217, 'Sb'  : 0.206, 'Te' : 0.206, 
                 'I'   : 0.198, 'Xe'  : 0.216, 'Cs'  : 0.167, 'Ba' : 0.149, 
                 'La'  : 0.200, 'Ce'  : 0.200, 'Pr'  : 0.200, 'Nd' : 0.200, 
                 'Pm'  : 0.200, 'Sm'  : 0.200, 'Eu'  : 0.200, 'Gd' : 0.200, 
                 'Tb'  : 0.200, 'Dy'  : 0.200, 'Ho'  : 0.200, 'Er' : 0.200, 
                 'Tm'  : 0.200, 'Yb'  : 0.200, 'Lu'  : 0.200, 'Hf' : 0.200, 
                 'Ta'  : 0.200, 'W'   : 0.200, 'Re'  : 0.200, 'Os' : 0.200, 
                 'Ir'  : 0.200, 'Pt'  : 0.175, 'Au'  : 0.166, 'Hg' : 0.155, 
                 'Tl'  : 0.196, 'Pb'  : 0.202, 'Bi'  : 0.207, 'Po' : 0.197, 
                 'At'  : 0.202, 'Rn'  : 0.220, 'Fr'  : 0.348, 'Ra' : 0.283, 
                 'Ac'  : 0.200, 'Th'  : 0.200, 'Pa'  : 0.200, 'U'  : 0.186, 
                 'Np'  : 0.200, 'Pu'  : 0.200, 'Am'  : 0.200, 'Cm' : 0.200, 
                 'Bk'  : 0.200, 'Cf'  : 0.200, 'Es'  : 0.200, 'Fm' : 0.200, 
                 'Md'  : 0.200, 'No'  : 0.200, 'Lr'  : 0.200, 'Rf' : 0.200, 
                 'Db'  : 0.200, 'Sg'  : 0.200, 'Bh'  : 0.200, 'Hs' : 0.200, 
                 'Mt'  : 0.200, 'Ds'  : 0.200, 'Rg'  : 0.200, 'Cn' : 0.200, 
                 'Uut' : 0.200, 'Fl'  : 0.200, 'Uup' : 0.200, 'Lv' : 0.200, 
                 'Uus' : 0.200, 'Uuo' : 0.200} 

import argparse
import math
import numpy as np
import itertools

def vdwfr_min_calc(alist1,alist2,crd1,crd2): #minimum vdw fraction
  natom1,natom2=len(alist1),len(alist2)
  vdwfr=np.empty((0),dtype=float)

  for i in range(natom1):
    for j in range(natom2):
      rij=crd1[i]-crd2[j]
      dist=np.sqrt(np.dot(rij,rij))
      sigma=_ATOMIC_RADII[alist1[i]]+_ATOMIC_RADII[alist2[j]]
      vdwfr1=dist/sigma
      vdwfr=np.append(vdwfr,vdwfr1)

  vdwfr_min=np.amin(vdwfr)
  return vdwfr_min

parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='calculate distance fraction with respect to vdw contact surface')
## args
parser.add_argument('-i', '--input', default='monomer.xyz', nargs='?', 
  help='input xyz file (following xyz format in jmol)')
#parser.add_argument('-probe_r', '--probe_radius', default=0.12, nargs='?', type=float,
#  help='target distance from VdW boundary of the monomer (nm)')
parser.add_argument('-n1', '--n1_xyz', default=1, nargs='?', type=int,
  help='number of atoms in monomer 1')
parser.add_argument('-n2', '--n2_xyz', default=1, nargs='?', type=int,
  help='number of atoms in monomer 2')
#parser.add_argument('-s', '--sigma_ratio', default=1.0, nargs='?', type=float,
#  help='the radius in unit of vdw sigma for rolling distance')
parser.add_argument('-o', '--output', default='monomer_h_', nargs='?', 
  help='output prefix for output xyz files')
parser.add_argument('args', nargs=argparse.REMAINDER)
#parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))

xyzfile = open(args.input,'r')
outfile = open(args.output,'w')
natom1=args.n1_xyz
natom2=args.n2_xyz
print("natoms1,natoms2 = ",natom1,natom2)
n1set=natom1+natom2+3 #number of lines of one dimer set in xyz file

#loop : should continue reading until EOF.
while True:
  next_n_lines=list(itertools.islice(xyzfile,n1set))
  if not next_n_lines: #EOF
    break
  else: #should exactly match number of lines for 1 dimer
    data1 = [line.split() for line in next_n_lines[2:natom1+2]]
    data2 = [line.split() for line in next_n_lines[natom1+2:n1set-1]]

    np_data1,np_data2=np.array(data1),np.array(data2)
    alist1,alist2=np_data1[:,0],np_data2[:,0]
    #convert unit: angstrom -> nm
    crd1,crd2=np.array(np_data1[:,1:4],np.float32)/10.0,np.array(np_data2[:,1:4],np.float32)/10.0

    vdwfr_min=vdwfr_min_calc(alist1,alist2,crd1,crd2)
    outfile.write('{:8.3f}\n'.format(vdwfr_min))

