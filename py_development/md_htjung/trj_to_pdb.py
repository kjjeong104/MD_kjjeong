# arg[1] : trajectory filename
# arg[2] : nth frame to get in [1:n_frames]
# arg[3] : topology file

import mdtraj as md
import sys
print(sys.argv)
t = md.load_frame(sys.argv[1],int(int(sys.argv[2])-1),sys.argv[3])
if len(t) == 0:
	raise ValueError(" wrong second argument!")
t.save_pdb("out.pdb")
