rm *log md* *.out

head -n -3 log > energy.log
sed -i '1,20d' energy.log
sed -i '/class/d' energy.log
awk ' NR % 2 == 0' energy.log > energy_pot.log
awk '{ print $1 }' energy_pot.log > energy_pot.out

import numpy as np
import mdtraj as md
import math 
#t = md.load('md_nvt.dcd', top='start_drudes.pdb') # gas
t = md.load('md_npt.dcd', top='start_drudes.pdb') # liq
vols = t.unitcell_lengths[:,0]**3
n_frames = len(vols)
np.mean(vols[int(n_frames/2):])
angles = t.top.select("name H9 or name C3 or name H10")
list_ang = np.reshape(angles, (-1,3))
list_ang[:,0], list_ang[:,1] = list_ang[:,1], list_ang[:,0].copy() 
result_ang = md.compute_angles(t, list_ang)
n_frames = len(result_ang)
np.mean(result_ang[int(n_frames/2):])*180.0/math.pi
np.std(result_ang[int(n_frames/2):])*180.0/math.pi

import numpy as np
pot = np.loadtxt('energy_pot.out')
n_frames = len(pot)
np.mean(pot[int(n_frames/2):])

choline
104.17

TFSI
280.15

TMPA
102.20


