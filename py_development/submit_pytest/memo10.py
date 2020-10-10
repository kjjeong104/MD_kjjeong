#!/usr/bin/python3.4

import sys
import numpy
import mdtraj as md

infile=open(sys.argv[1],'r')
outfile=open(sys.argv[2],'w')

for line in infile:
  outfile.write(line)
  outfile.write("hello world\n")

infile.close()
outfile.close()

