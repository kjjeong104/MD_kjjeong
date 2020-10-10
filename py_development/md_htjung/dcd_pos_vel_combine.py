import MDAnalysis as mda
import copy
import sys

fpos=sys.argv[1]
fvel=sys.argv[2]
fpdb=sys.argv[3]
u = mda.Universe(fpdb,fpos)
v = mda.Universe(fpdb,fvel)

# set the has_velocities flag on the Timestep
u.trajectory.ts.has_velocities = True

# assign velocities to u
real = u.select_atoms("all and not name D*")
ww = mda.Writer("traj.trr", len(real))
print(" total {} frames".format(len(u.trajectory)))
for iframe in range(len(u.trajectory)):
	#print(iframe)
	ts=u.trajectory[iframe]
	tsv=v.trajectory[iframe]
	u.atoms.velocities = copy.copy(v.atoms.positions)
	#print(v.atoms.positions[0]) #unit = A/ps for mda, nm/ps for mdtraj
	if iframe == 0:
		#print(v.atoms.positions[0])
		real.write("conf.gro") # unit = nm, nm/ps
	ww.write(real)
ww.close()
