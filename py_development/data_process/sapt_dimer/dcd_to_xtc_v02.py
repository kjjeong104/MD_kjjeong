#!/usr/bin/env python
# program to convert dcd trajectroy to xtc file
#inputs: pdbfile, traj dcd file.
#output: target xtc file

import math
import sys
import numpy
#import mdtraj as md
import MDAnalysis as mda
import copy

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  outfile = sys.argv[3]

  #input 1 : load traj. (big file)
  u=mda.Universe(topfile,trjfile)
  real = u.select_atoms("all")
  ww=mda.Writer(outfile, len(real))
  print(" total {} frames".format(len(u.trajectory)))

  for iframe in range(len(u.trajectory)):
    ts=u.trajectory[iframe]
    if iframe == 0:
      real.write(topfile+".gro")
    ww.write(real)
  ww.close()
  #traj=md.load(trjfile,top=topfile)
  #traj.save_xtc(outfile)

if __name__ == "__main__": main()

