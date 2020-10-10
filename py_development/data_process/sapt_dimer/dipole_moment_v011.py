#!/usr/bin/env python
#program to read trajectory file, net charge value table and evaluate dipole moment 
#inputs: pdbfile, traj (dcd or xtc) file, charge table file
#quick calculator : urea bl
#v02 : can process trjaectory

#import math
import sys
import numpy
import mdtraj as md
import argparse

#def find_nearest_index(array, value):
    #array = np.asarray(array)
#    idx = (numpy.abs(array - value)).argmin()
#    return idx
conv_D_enm = 0.020819434
conv_deg_rad= numpy.pi/180.0
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
parser.add_argument('-u', '--urefile', default='urea.txt', nargs='?', 
  help='optional urea overwritten geometric variables')
parser.add_argument('-o', '--outfile', default='dipole.dat', nargs='?', 
  help='dipole moment history for the trajectory')
args = parser.parse_args()

#def urea_constructor(topfile,urefile):
#  newtraj=md.load(topfile)
#
#  newtraj.xyz[0][cindex]=[0.0,0.0,0.0] #draw C atom
#  newtraj.xyz[0][oindex]=[0.0,blco,0.0] #draw O atom
#  newtraj.xyz[0][nindex[0]],newtraj.xyz[0][nindex[1]]= #draw N atoms
   
#  return newtraj

#main fxn
def main():
  #Part to load coordinate file
  topfile = args.topfile
  trjfile = args.trjfile
  chgfile = args.chgfile
  #urefile = args.urefile
  outfile = open(args.outfile,'w')
  carray=numpy.loadtxt(chgfile)
 
  #input 1 : load surf traj. (big file)
  #t=md.load(topfile)
  traj=md.load(trjfile,top=topfile)
  nstep=traj.n_frames
  #moments=md.dipole_moments(traj,carray)
  moments=md.dipole_moments(traj,carray)
  #print(moments)
  for i in range(nstep):
    outfile.write("{:8.3f} D\n".format(numpy.linalg.norm(moments[i])/conv_D_enm))
  #print(str(numpy.linalg.norm(moments[0])/conv_D_enm)+' D')

  #optional : if urea overwritten geometric variables exist, parse that.
 # newtraj=urea_constructor(topfile,urefile)
 # moments=md.dipole_moments(newtraj,carray)
 # print(moments)
 # print(str(numpy.linalg.norm(moments[0])/conv_D_enm)+' D')

  #numpy.savetxt('boxdim'+trjfile[:-4]+'.dat',boxinfo)
  #idx=find_nearest_index(boxinfo,boxdemand)
  #snapshot=traj[idx]
  #topology=snapshot.topology
  #drude particle exclusion
  #snapshot=snapshot.atom_slice([atom.index for atom in topology.atoms if ('D' not in atom.name) ])
  #snapshot.save_pdb(outfile)

if __name__ == "__main__": main()

