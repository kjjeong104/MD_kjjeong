#!/usr/bin/env python
#inputs: fragmental pdb file
#output: new fragmental pdb file

import sys
import math
#main fxn
def main():
  #Part to load coordinate file
  inpfile=open(sys.argv[1],'r')
  outfile=open(sys.argv[2],'w')
  inc=int(input("increment of 2nd column? \n"))
  period=int(input("Period of residues? \n"))
  comm=input("Want 'shift' mode or 'reset' mode for atom index?\n")
  
  if comm=='shift':
    for line in inpfile:
      num=int(line[22:26]) #col 22-25 of pdb file format
      num2=num+inc
      #num3=math.ceil(num2)
      outfile.write(line[0:22]+'{:4}'.format(num2)+line[26:])
  elif comm=='reset':
    num2=inc
    for line in inpfile:
      num2+=(1/period)
      num3=math.ceil(num2)
      outfile.write(line[0:22]+'{:4}'.format(num3)+line[26:])

if __name__ == "__main__": main()

