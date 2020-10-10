#!/home/kjeong23/softwares/bin/python3.4
#program to pick a snapshot with suitable box size
#from npt simulation trajectory
# program to convert dcd trajectroy to xtc file
#inputs: pdbfile, traj dcd file.
#output: target xtc file

#import math
import sys
import numpy
import mdtraj as md

def find_nearest_index(array, value):
    #array = np.asarray(array)
    idx = (numpy.abs(array - value)).argmin()
    return idx

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  outfile = sys.argv[3]
  boxdemand=float(input("How much is the box size demanded(in nm)? ex) 4.34 \n"))
  teq=int(input("How many initial frames do you want to cut as equilibration? ex) 5000 \n"))


  #input 1 : load surf traj. (big file)
  traj=md.load(trjfile,top=topfile)
  traj=traj[teq:]
  nstep=traj.n_frames
  boxinfo=traj.unitcell_lengths[:,0] #assuming cubic box.
  #traj.save_xtc(outfile)

  idx=find_nearest_index(boxinfo,boxdemand)
  snapshot=traj[idx]
  topology=snapshot.topology
  #drude particle exclusion
  snapshot=snapshot.atom_slice([atom.index for atom in topology.atoms if ('D' not in atom.name) ])

  snapshot.save_pdb(outfile)

if __name__ == "__main__": main()

