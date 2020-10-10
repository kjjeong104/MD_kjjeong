#!/usr/bin/env python
#inputs: fragmental pdb file
#output: new fragmental pdb file

import sys

#main fxn
def main():
  #Part to load coordinate file
  inpfile=open(sys.argv[1],'r')
  outfile=open(sys.argv[2],'w')
  inc=int(input("increment of 2nd column? \n"))
  comm=input("Want 'shift' mode or 'reset' mode for atom index?\n")
  
  if comm=='shift':
    for line in inpfile:
      num=int(line[6:11]) #col 7-11 of pdb file format
      num2=num+inc
      outfile.write(line[0:6]+'{:5}'.format(num2)+line[11:])
  elif comm=='reset':
    num2=inc
    for line in inpfile:
      num2+=1
      outfile.write(line[0:6]+'{:5}'.format(num2)+line[11:])

if __name__ == "__main__": main()

