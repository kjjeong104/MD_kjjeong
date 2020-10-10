#!/usr/bin/env python
# program to convert dcd trajectroy to xtc file
#inputs: pdbfile, traj dcd file.
#output: target xtc file

import math
import sys
import numpy
import mdtraj as md

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  outfile = sys.argv[3]

  tskip=int(input("Once in how many frames do you want to take? ex) 10 \n"))
  teq=int(input("How many initial frames do you want to cut as equilibration?(before counting tskip) ex) 2000 \n"))
  #input 1 : load surf traj. (big file)
  traj=md.load(trjfile,top=topfile)
  print("#frames original : ",traj.n_frames)
  traj=traj[teq:]
  traj=traj[::tskip]
  print("#frames after trimming : ",traj.n_frames)
  traj.save_dcd(outfile)

if __name__ == "__main__": main()

