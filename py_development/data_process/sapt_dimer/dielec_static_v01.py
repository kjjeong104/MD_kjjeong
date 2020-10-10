#!/usr/bin/env python
#program to read trajectory file, net charge value table 
#and evaluate static dielectric constant.
#inputs: pdbfile, traj (dcd or xtc) file, charge table file
#charge table file should have format that list 1 molecule name,
# then charge values of 1 molecule atoms in the same order of pdb atom order.

#import math
import sys
import numpy
import mdtraj as md
import argparse

#conv_D_enm = 0.020819434
#conv_deg_rad= numpy.pi/180.0
parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='dipole moment calculator')
# args
parser.add_argument('-s', '--topfile', default='', nargs='?', 
  help='topology pdb file')
parser.add_argument('-f', '--trjfile', default='', nargs='?', 
  help='trajectory file')
parser.add_argument('-c', '--chgfile', default='carray.txt', nargs='?', 
  help='charge value array file')
parser.add_argument('-temp', '--temperature', default=300.0, nargs='?', 
  help='temperature of the run in Kelvin')
parser.add_argument('-o', '--outfile', default='static_dielectric.dat', nargs='?', 
  help='dipole moment history for the trajectory')
parser.add_argument('-b', '--begin',default=0,nargs='?',
  help='initial cut frames for equilibration ')
args = parser.parse_args()

#read topology file and 1 molecule chgfile, to construct full charge-array
def construct_carray(topfile,chgfile): 
  t=md.load(topfile)
  top=t.topology
  #use python dictionary attribute to store monomer point charge information, and append to the grand charge-array
  #read chgfile
  f=open(chgfile,'r')
  f1=f.readlines()
  resmat=[i for i, x in enumerate(f1) if x[0].isalpha()==True] #find out linenumbers in chgfile with resnames.
  cdict={}
  for i in range(len(resmat)):
    x=resmat[i]
    rname=f1[x].rstrip()
    if i+1<len(resmat):
      y=resmat[i+1]
      onemc=numpy.array([float(z) for z in f1[x+1:y]])
    else:
      onemc=numpy.array([float(z) for z in f1[x+1:]])
    cdict.update({rname : onemc })
    print(' monomer detected. name ',rname,' charge ',numpy.sum(onemc))
  #transcript the topology file residue info to the grand charge-array
  #print(cdict)
  carray=numpy.empty(0)
  for i in range(t.n_residues):
    rname=top.residue(i).name
    carray=numpy.append(carray,cdict[rname])
  return carray

#main fxn
def main():
  #Part to load coordinate file
  topfile = args.topfile
  trjfile = args.trjfile
  chgfile = args.chgfile
  #urefile = args.urefile
  outfile = open(args.outfile,'w')
  #carray=numpy.loadtxt(chgfile)
  carray=construct_carray(topfile,chgfile)
  T,teq = float(args.temperature),int(args.begin)
  #input 1 : load surf traj. (big file)
  #t=md.load(topfile)
  traj=md.load(trjfile,top=topfile)
  traj=traj[teq:]
  nstep=traj.n_frames
  print(nstep, ' frames ')
  #moments=md.dipole_moments(traj,carray)
  epsil=md.static_dielectric(traj,carray,T)
  outfile.write("Dielectric constant {:8.3f}\n".format(epsil))
  #print(epsil)
  #print(moments)
  #for i in range(nstep):
  #  outfile.write("{:8.3f} D\n".format(numpy.linalg.norm(moments[i])/conv_D_enm))
  #optional : if urea overwritten geometric variables exist, parse that.
 # newtraj=urea_constructor(topfile,urefile)
 # moments=md.dipole_moments(newtraj,carray)
 # print(moments)
 # print(str(numpy.linalg.norm(moments[0])/conv_D_enm)+' D')
  outfile.close()

if __name__ == "__main__": main()

