import numpy as np
import mdtraj as md
import math 
import sys

phase=sys.argv[1]
mass=sys.argv[2]

if 'gas' in phase:
	t = md.load('md_nvt.dcd', top='start_drudes.pdb') # gas
elif 'liq' in phase:
	t = md.load('md_npt.dcd', top='start_drudes.pdb') # liq
else:
	raise ValueError("correct phase? {}".format(str(phase)))

# calculate average density
vols = t.unitcell_lengths[:,0]**3
n_frames = len(vols)
print(" total {} frames".format(n_frames))
avg_vol = np.mean(vols[int(n_frames/2):])
np.savetxt("density.log",np.column_stack((float(mass)*10.0/6.022/vols,t.unitcell_lengths[:,0])))
density = float(mass)*10.0/6.022/avg_vol
print(" avg. density = {}".format(density))
# find a pdb file close to average density
idx_pdb = np.argmin(np.abs(vols - avg_vol))
print(" pdb index for average density = {} ".format(idx_pdb))
print(" the pdb file has {}".format(t.unitcell_lengths[idx_pdb,0]))

angles = t.top.select("name H9 or name C3 or name H10")
list_ang = np.reshape(angles, (-1,3))
list_ang[:,0], list_ang[:,1] = list_ang[:,1], list_ang[:,0].copy() 
result_ang = md.compute_angles(t, list_ang)
n_frames = len(result_ang)
ang_avg = np.mean(result_ang[int(n_frames/2):])*180.0/math.pi
ang_std = np.std(result_ang[int(n_frames/2):])*180.0/math.pi
print(" H-C(N)-H angle = {} +- {}".format(ang_avg,ang_std))

pot = np.loadtxt('energy_pot.out')
n_frames = len(pot)
print(n_frames/2)
avg_pot = np.mean(pot[int(n_frames/2)+1:])
print(" avg. potential energy = {}".format(avg_pot))


