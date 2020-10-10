#!/home/kjeong23/softwares/bin/python3.4
# program to collect and convert psi4 SAPT calculation output files to .bohr file compatible to fitting code
#inputs: .xyz files for atom names, psi4 .out files, max number of output index
#output: .bohr file

import math
import sys
import numpy
import os
import os.path

#Parameters
angtobohr=1.889725989 #conversion factor angstrom -> bohr
en_factor=1000.0 #Ha->mH
disp_scale=1.06 #dispersion energy scaling factor
dftmethod='pbe0'
reverse=True #swap the monomer 1 and monomer 2
#dictionary of necessary energy components
ecompstr=['E1pol', 'E1exch', 'E1exch(S2)', 'E2ind(unc)', 'E2ind', \
'E2ind-exch', 'E2disp(unc)', 'E2disp', 'E2disp-exch(unc)', 'E2disp-exch', \
'E1tot', 'E2tot', 'E1tot+E2tot', 'E2ind[B->A]', 'E2ind[A->B]', 'E2exchind_BA', 'E2exchind_AB', 'dhf']
psi4ecomp={'Elst1,r   ':0, 'Exch1   ':1, 'Exch1(S^2)   ':2, }


def xyz_angtobohr(line)
  lsplit=line.split()
  xb,yb,zb=float(lsplit[1]),float(lsplit[2]),float(lsplit[3])
  xb,yb,zb="{:15.8f}".format(xb*angtobohr),"{:15.8f}".format(yb*angtobohr),"{:15.8f}".format(zb*angtobohr)
  addrow=' '+xb+' '+yb+' '+zb
  return addrow

#function to write 1 paragraph of bohr file, reading 1 psi4 output
def psi4_bohr_paragraph(psifile,aname1,aname2):
  paragraph=[] #empty string list
  natom1,natom2=len(aname1),len(aname2)
  ecomp=numpy.zeros(len(ecompstr))

  while True:
    line=psifile.readline()
    if line=='':#EOF
      break
    if 'molecule dimer ' in line: #coord block
      line=psifile.readline() #skip 1 line
      #monomer 1 coord reading and writing
      paragraph.append(natom1)
      for i in range(natom1):
        row=aname1[i]
        line=psifile.readline()
        row+=xyz_angtobohr(line)
        paragraph.append(row)
      #monomer 2 coord reading and writing
      line=psifile.readline() #skip 2 lines
      line=psifile.readline()
      paragraph.append(natom2)
      for i in range(natom2):
        row=aname2[i]
        line=psifile.readline()
        row+=xyz_angtobohr(line)
        paragraph.append(row)

    #energy component reading. Be careful of errorneous values
    elif 'Elst1,r   '  in line: ecomp[0]=float(line.split()[1])
    elif 'Exch1   '  in line: ecomp[1]=float(line.split()[1])
    elif 'Exch1(S^2)   '  in line: ecomp[2]=float(line.split()[1])
    elif in line: ecomp[3]=float(line.split()[1])
    elif in line: ecomp[4]=float(line.split()[1])
in line: ecomp[5]=float(line.split()[1])
in line: ecomp[6]=float(line.split()[1])
in line: ecomp[7]=float(line.split()[1])
in line: ecomp[8]=float(line.split()[1])

  #after finishing energy component reading, write paragraph
  for i in range(len(ecompstr)):
    row=ecompstr[i]+' '+ecomp[i]
    paragraph.append(row)
  #1 empty line
  paragraph.append(' ')

  return paragraph #string list

#main fxn
def main():
  #Part to load and put information
  mon12=input("Name of two monomers? ex) EMIM Cl\n")
  msplit=mon12.split()
  namemon1,namemon2=msplit[0],msplit[1]
  nfiles=int(input("Number of input files used in this data set? ex) 1000\n"))
  an1file=namemon1+'.xyz'
  an2file=namemon2+'.xyz'
  bohrfile = sys.argv[1]

  #get list of anames
  aname1,aname2=[],[]
  for line in an1file:
    lsplit=line.split()
    if reverse==True: #swap the monomer 1 and monomer 2
      aname2.append(lsplit[0])
    else:
      aname1.append(lsplit[0])
  for line in an2file:
    lsplit=line.split()
    if reverse==True: #swap the monomer 1 and monomer 2
      aname1.append(lsplit[0])
    else:
      aname2.append(lsplit[0])

  for i in range(nfiles): #loop for reading each psi4 output files to write dataset in .bohr file
    psifile=str(namemon1+'_'+namemon2+'_'+dftmethod+'_'+str(i)+'.out')
    pathtest='./'+psifile
    if os.path.isfile(pathtest) and os.access(pathtest, os.R_OK): #output file exists
      paragraph=psi4_bohr_paragraph(psifile,aname1,aname2)
      for row in paragraph:
        bohrfile.write(row)
      print('Psi4 output file index {} interpreted'.format(i))

  bohrfile.close()

if __name__ == "__main__": main()

