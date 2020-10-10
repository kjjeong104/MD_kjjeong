#!/home/kjeong23/softwares/bin/python3.4
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

  #input 1 : load surf traj. (big file)
  traj=md.load(trjfile,top=topfile)
  traj.save_xtc(outfile)

if __name__ == "__main__": main()

