import MDAnalysis as mda
import copy
import sys

fpos=sys.argv[1]
fpdb=sys.argv[2]
u = mda.Universe(fpdb,fpos)

# set the has_velocities flag on the Timestep
#u.trajectory.ts.has_velocities = True

# assign velocities to u
real = u.select_atoms("all and not name D*")
#real = u.select_atoms("all")
ww = mda.Writer("traj_noD.trr", len(real))
print(" total {} frames".format(len(u.trajectory)))
for iframe in range(len(u.trajectory)):
	#print(iframe)
	ts=u.trajectory[iframe]
	#print(v.atoms.positions[0]) #unit = A/ps
	if iframe == 0:
		#print(v.atoms.positions[0])
		real.write("conf.gro") # unit = nm, nm/ps
	ww.write(real)
ww.close()
