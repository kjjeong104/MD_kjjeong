#!/usr/bin/env python3
# ver 0.1 - coding python by Hyuntae Jung on 10/10/2018

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
                 'Uus' : 0.200, 'Uuo' : 0.200} #in nm 

import argparse
parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='generate grid points for electric field')
## args
parser.add_argument('-i', '--input', default='monomer.xyz', nargs='?', 
  help='input xyz file (following xyz format in vmd)')
parser.add_argument('-probe_r', '--probe_radius', default=0.12, nargs='?', type=float,
  help='target distance from VdW boundary of the monomer (nm)')
parser.add_argument('-d','--dist_points',default=0.05, nargs='?', type=float,
  help='threshold distance between probe points (lower bound, nm)')
parser.add_argument('-n', '--n_xyz', default=1000, nargs='?', type=int,
  help='the number of points for output')
parser.add_argument('-sigma', '--sigma_ratio', default=1.0, nargs='?', type=float,
  help='the radius in unit of vdw sigma for rolling distance')
parser.add_argument('-o', '--output', default='monomer_h_', nargs='?', 
  help='output prefix for output xyz files')
parser.add_argument('args', nargs=argparse.REMAINDER)
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
## read args
args = parser.parse_args()
## Check arguments for log
print(" input arguments: {0}".format(args))


#```
#Compute the coordinates of points on a sphere using the Golden Section Spiral algorithm.
#n_pts : Number of points to generate on the sphere
#```

import math
import numpy as np

f = open(args.input,'r')
lines = f.readlines()
data = [line.split() for line in lines[2:]]
dim1 = len(data) 
atom_mapping = np.arange(dim1, dtype=np.int32)
n_atoms = dim1

np_data = np.array(data)
xyzlist = np.array(np_data[:,1:4],np.float32)
atoms = np_data[:,0]
sigma=args.sigma_ratio
atom_radii = [_ATOMIC_RADII[atom] for atom in atoms]
atom_radii = (np.array(atom_radii, np.float32) + args.probe_radius)*10.0*sigma #conversion to angstroms

def generate_sphere_points(n_points): 
  inc = math.pi * (3.0 - np.sqrt(5.0))
  offset = 2.0 / float(n_points)
  sphere_points = np.empty((n_points,3),dtype=np.float32)

  for i in range(n_points):
    y = i*offset - 1.0 + (offset / 2.0)
    r = np.sqrt(1.0 - y*y)
    phi = i * inc
    sphere_points[i] = np.array([math.cos(phi) * r, y, math.sin(phi) * r], dtype=np.float32)

  return sphere_points

n_total_points = args.n_xyz
n_sphere_points = int(n_total_points/n_atoms*3.0) #doesn't need too many candidate points
sphere_points = generate_sphere_points(n_sphere_points)

#rolling surface construction algorithm
#step 1 : generate centered sphere points for all atoms. during that generation,
#test if too closer to other atoms. if closer to other atoms, don't accept it
#step 2: among the chosen probe points, sort out points according to probe-probe distance threshold
#only leave designated number of points. we have to evenly choose points around the space
total_probes1=np.empty((0,3),dtype=float)
probe_parent=[]
for i in range(n_atoms):
  atom_radius_i = atom_radii[i]
  sphere_cand1=xyzlist[i] + atom_radius_i*sphere_points
  for point in sphere_cand1:
    register=True
    for j in range(n_atoms):
      if (i==j):
        continue
      else:
        r_j = xyzlist[j]
        r_pj= point - r_j
        r_pj2=np.dot(r_pj,r_pj)
        if(r_pj2<atom_radii[j]):
          register=False
          break
    if(register):
      total_probes1=np.vstack((total_probes1,point))
      probe_parent.append(i)

for point in total_probes1:
  print('X {:8.3f} {:8.3f} {:8.3f}'.format(point[0],point[1],point[2]))

