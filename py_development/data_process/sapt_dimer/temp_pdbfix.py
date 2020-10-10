#!/usr/bin/env python
#inputs: pdb file
#output: pdb file
#correct strange spacing of pdb file (FPMD traj) into readable position

import sys

#main fxn
def main():
  #Part to load coordinate file
  inpfile=open(sys.argv[1],'r')
  outfile=open(sys.argv[2],'w')
  #inc=int(input("increment of 2nd column? \n"))

  #correct spacing:
  '''
HETATM    1  O2a ure A   1      28.649  14.336  38.706  1.00  0.00           O 
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom. 
  '''
  for line in inpfile:
    #6 5 5 5 0 5 : 6 11 16 21 21 26
    #8 8 8 6 6 2 : 34 42 50 56 62 64 : additional empty space demand: 4 0 0 0 0 10
    # 1~6, 8~12, 13~17,18~22,empty,24~28
    # 30~37, 39~46, 48~55, 57~62, 64~69,77~78 
    seg1,seg2,seg3,seg4,seg5,seg6=line[0:6],line[7:12],line[12:17],line[17:22],line[22:22],line[23:28]
    seg7,seg8,seg9,seg10,seg11=line[29:37],line[38:46],line[47:55],line[56:62],line[63:69]
    seg12=line[76:]
    #if len(line)<78:
    #  seg12=line[76:77]+'  '
    #else:
    #  seg12=line[76:78]+' '
    outfile.write(seg1+seg2+seg3+seg4+seg5+seg6+'    '+seg7+seg8+seg9+seg10+seg11+'     '+'     '+seg12)

if __name__ == "__main__": main()

