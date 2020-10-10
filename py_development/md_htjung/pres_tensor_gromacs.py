## reding trajectory to calculate pressure tensor in GROMACS

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

u = md.Universe('conf2.gro','traj.trr')
n_frames = len(u.trajectory)
mass_array = u.atoms.masses # unit = [dalton] in GROMACS
#print(" mass array = {} [kg/mol]".format(mass_array))

#for ts in u.trajectory:
for i in range(n_frames):
	ts=u.trajectory[i]
	atomic_vel1 = np.array(ts._velocities) # unit = [A/ps] 
	if i+1 == n_frames:
		atomic_vel = atomic_vel1
	else:
		ts=u.trajectory[i+1]
		atomic_vel2 = np.array(ts._velocities) # unit = [A/ps] 
		# to reduce energy error came from integrator
		#  we use average velocity, which is the similar way with GROMACS.
		atomic_vel = 0.5*(atomic_vel1+atomic_vel2)
	
	ts=u.trajectory[i]
	atomic_pos = np.array(ts._pos) # unit = [A]
	atomic_for = np.array(ts._forces) # unit = [kj/(A*mole)]
	# calculate kinetic energy term for P_xy
	kin_e = np.empty(6) 
	kin_e[0] = 0.5*np.sum(mass_array*(atomic_vel[:,0]*atomic_vel[:,0]))/100 # xx
	kin_e[1] = 0.5*np.sum(mass_array*(atomic_vel[:,0]*atomic_vel[:,1]))/100 # xy
	kin_e[2] = 0.5*np.sum(mass_array*(atomic_vel[:,0]*atomic_vel[:,2]))/100 # xz
	kin_e[3] = 0.5*np.sum(mass_array*(atomic_vel[:,1]*atomic_vel[:,1]))/100 # yy
	kin_e[4] = 0.5*np.sum(mass_array*(atomic_vel[:,1]*atomic_vel[:,2]))/100 # yz
	kin_e[5] = 0.5*np.sum(mass_array*(atomic_vel[:,2]*atomic_vel[:,2]))/100 # zz
	## Be aware that the calculation result would be different with gmx energy
	#  1. integrator (in case of Verlet, pos and vel are obtained or stored at different time step in trajectory file)
	#     the positions are a half time step later than the velocities, so potential and kinetic energies are not at the same timestep 
	#     See details: https://github.com/pandegroup/openmm/issues/1474
	#  2. epecially, single precision problem 
	#     (reading trajectory in mdtraj is not perfectly same as how gromacs module works, thus conversion error occurs < 1.0%)
	#   ex) -1.47786 (original) --> 14.778567 (reading)
	#   ex) -33.3855 (original) --> 3.338553  (reading)
	total_kin_e = kin_e[0]+kin_e[3]+kin_e[5] 
	print(" kin. energy by mannual= {} kJ/mol".format(total_kin_e))
	#print(" kin. energy by OpenMM = {}".format(str(state.getKineticEnergy()))) 
	
	## check center of mass velocity
	reduce_atomic_vel = np.empty(np.shape(atomic_vel))
	n_atoms = len(mass_array)
	for i_atom in range(n_atoms):
		reduce_atomic_vel[i_atom] = mass_array[i_atom] * atomic_vel[i_atom]
	com_vel = np.sum(reduce_atomic_vel, axis=0)/np.sum(mass_array)
	com_kin_e = 0.5*np.sum(mass_array)*np.sum(com_vel*com_vel)
	if com_kin_e*100/total_kin_e > 1.0:
		print(" COM kinetic energy = {} (should be removed!)".format(com_kin_e))
		print(total_kin_e - com_kin_e)
	
	## calculate virial term
	box_3dim = ts.dimensions[:3]
	pair_dis = pairwise_row_diff_slicing(atomic_pos) # difference between any two vectors in atomic_pos without duplicates
	pair_dfor = pairwise_row_diff_slicing(atomic_for) # difference between any two vectors in atomic_pos without duplicates
	pair_dis_short = pair_dis - np.around(pair_dis/box_3dim)*box_3dim # reduce into the shortest vectors due to periodic boundary condition
	
	virial = np.empty(6) 
	virial[0] = -np.sum(pair_dis_short[:,0]*pair_dfor[:,0])/2. # xx
	virial[1] = -np.sum(pair_dis_short[:,0]*pair_dfor[:,1])/2. # xy
	virial[2] = -np.sum(pair_dis_short[:,0]*pair_dfor[:,2])/2. # xz
	virial[3] = -np.sum(pair_dis_short[:,1]*pair_dfor[:,1])/2. # yy
	virial[4] = -np.sum(pair_dis_short[:,1]*pair_dfor[:,2])/2. # yz
	virial[5] = -np.sum(pair_dis_short[:,2]*pair_dfor[:,2])/2. # zz
	#print(" virial = {} kJ/mol".format(virial)) # unit = [kJ/mol]
	box_vol = np.product(box_3dim)

	pressure_tensor = 2.*(kin_e - virial)*100000/(box_vol*6.0221417930) # unit: bar
	#print(" pressure tensor = {} bar".format(pressure_tensor))
	print(" intantaneous pressure (iso) = {} bar".format(
			(pressure_tensor[0]+pressure_tensor[3]+pressure_tensor[5])/3.)) # finally, %5 error with GROMACS gmx energy result
	
#box_vol = state.getPeriodicBoxVolume()
#box_vol = box_vol.in_units_of(meter**3)._value
#print(virial)
#pressure_tensor = (kin_e + virial)/(box_vol*6.0221417930e23)/1000 # unit: MPa
#print("est. pressure = {}".format((pressure_tensor[0]+pressure_tensor[3]+pressure_tensor[5])/3))
#print(pressure_tensor)
#print(str(state.getPotentialEnergy()))
#for j in range(system.getNumForces()):

