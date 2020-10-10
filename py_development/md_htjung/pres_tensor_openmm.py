
import MDAnalysis as md
import numpy as np

def pairwise_row_diff_slicing(a):
    n = len(a)
    N = n*(n-1)//2
    idx = np.concatenate(( [0], np.arange(n-1,0,-1).cumsum() ))
    start, stop = idx[:-1], idx[1:]
    out = np.empty((N,a.shape[1]),dtype=a.dtype)
    for j,i in enumerate(range(n-1)):
        out[start[j]:stop[j]] = a[i+1:] - a[i,None]
    return out

u = md.Universe('conf.gro','traj.trr')
n_frames = len(u.trajectory)
mass_array = u.atoms.masses # unit = [g/mol]
#for ts in u.trajectory:
ts=u.trajectory[1]
atomic_pos = np.array(ts._pos) # unit = [A]
atomic_vel = np.array(ts._velocities) # unit = [A/ps] = 1/10 * [nm/ps]
atomic_for = np.array(ts._forces) # unit = [kj/(nm*mole)]

# calculate kinetic energy term for P_xy
kin_e = np.empty(6) 
# 1000 is for unit conversion of kilojoule
# if you use nm/ps for velocity, remove *10 factor
kin_e[0] = 0.5*np.sum(mass_array*(atomic_vel[:,0])**2)*10/1000 # xx
kin_e[1] = 0.5*np.sum(mass_array*(atomic_vel[:,0]*atomic_vel[:,1]))*10/1000 # xy
kin_e[2] = 0.5*np.sum(mass_array*(atomic_vel[:,0]*atomic_vel[:,2]))*10/1000 # xz
kin_e[3] = 0.5*np.sum(mass_array*(atomic_vel[:,1]*atomic_vel[:,1]))*10/1000 # yy
kin_e[4] = 0.5*np.sum(mass_array*(atomic_vel[:,1]*atomic_vel[:,2]))*10/1000 # yz
kin_e[5] = 0.5*np.sum(mass_array*(atomic_vel[:,2]*atomic_vel[:,2]))*10/1000 # zz
print(kin_e)
# Be aware: the small difference occurs due to Verlet integrator (e.g. )
#  the positions are a half time step later than the velocities, so potential and kinetic energies are not at the same timestep 
# See details: https://github.com/pandegroup/openmm/issues/1474
#print(" kin. energy by mannual= {}".format((kin_e[0]+kin_e[3]+kin_e[5]))) 
#print(" kin. energy by OpenMM = {}".format(str(state.getKineticEnergy()))) 

# calculate virial term
box_3dim = ts.dimensions[:3]
pair_dis = pairwise_row_diff_slicing(atomic_pos) # difference between any two vectors in atomic_pos without duplicates
pair_dfor = pairwise_row_diff_slicing(atomic_for) # difference between any two vectors in atomic_pos without duplicates
#pair_dis_short = pair_dis - np.around(pair_dis/box_3dim)*box_3dim # reduce into the shortest vectors due to periodic boundary condition
virial = np.empty(6) 
virial[0] = -np.sum(atomic_pos[:,0]*atomic_for[:,0])/2 # xx
virial[1] = -np.sum(atomic_pos[:,0]*atomic_for[:,1])/2 # xy
virial[2] = -np.sum(atomic_pos[:,0]*atomic_for[:,2])/2 # xz
virial[3] = -np.sum(atomic_pos[:,1]*atomic_for[:,1])/2 # yy
virial[4] = -np.sum(atomic_pos[:,1]*atomic_for[:,2])/2 # yz
virial[5] = -np.sum(atomic_pos[:,2]*atomic_for[:,2])/2 # zz

box_vol = state.getPeriodicBoxVolume()
box_vol = box_vol.in_units_of(meter**3)._value
print(virial)
pressure_tensor = (kin_e + virial)/(box_vol*6.0221417930e23)/1000 # unit: MPa
print("est. pressure = {}".format((pressure_tensor[0]+pressure_tensor[3]+pressure_tensor[5])/3))
print(pressure_tensor)
#print(str(state.getPotentialEnergy()))
#for j in range(system.getNumForces()):

