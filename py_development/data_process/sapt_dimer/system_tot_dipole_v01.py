#!/usr/bin/env python
#program to read trajectory file, net charge value table and evaluate dipole moment 
#inputs: pdbfile, traj (dcd or xtc) file, charge table file
#v02 : can process trajectory : get distribution, ave, fluctuation.
#sysv01 version : calculate system total dipole moment.

#import math
import sys
import numpy
import mdtraj as md
import argparse

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
parser.add_argument('-o', '--outfile', default='sys_dipole', nargs='?', 
  help='dipole moment histogram for the trajectory')
parser.add_argument('-b', '--begin',default=0,nargs='?',
  help='initial cut frames for equilibration ')
args = parser.parse_args()
#read topology file and 1 molecule chgfile, to construct full charge-array
#also construct dictionary of residue indices.
def construct_carray_midict(topfile,chgfile): 
  t=md.load(topfile)
  top=t.topology
  #use python dictionary attribute to store monomer point charge information, and append to the grand charge-array
  #read chgfile
  f=open(chgfile,'r')
  f1=f.readlines()
  resmat=[i for i, x in enumerate(f1) if x[0].isalpha()==True] #find out linenumbers in chgfile with resnames.
  cdict,midict={},{}
  for i in range(len(resmat)):
    x=resmat[i]
    rname=f1[x].rstrip()
    if i+1<len(resmat):
      y=resmat[i+1]
      onemc=numpy.array([float(z) for z in f1[x+1:y]])
    else:
      onemc=numpy.array([float(z) for z in f1[x+1:]])
    cdict.update({rname : onemc })
    midict.update({rname : numpy.empty(0)})
    print(' monomer detected. name ',rname,' charge ',numpy.sum(onemc))
  #transcript the topology file residue info to the grand charge-array
  #print(cdict)
  carray=numpy.empty(0)
  for i in range(t.n_residues):
    rname=top.residue(i).name
    midict[rname]=numpy.append(midict[rname],i)
    carray=numpy.append(carray,cdict[rname])
  return carray,cdict,midict

#main fxn
def main():
  #Part to load coordinate file
  topfile = args.topfile
  trjfile = args.trjfile
  chgfile = args.chgfile
  #urefile = args.urefile
  outstr = args.outfile
  carray,cdict,midict=construct_carray_midict(topfile,chgfile)
  teq = int(args.begin)
  #input 1 : load surf traj. (big file)
  #t=md.load(topfile)
  traj=md.load(trjfile,top=topfile)
  topology = traj.topology
  traj=traj[teq:]
  nstep=traj.n_frames
  print(nstep, ' frames ')
  
  ##calculate individual molecular dipole moments.
  ##For each species, collect differently.
  #calculate system total dipole moment.
  outfile=outstr+'.dat'
  moments=md.dipole_moments(traj,carray)
  moments_normD=numpy.linalg.norm(moments,axis=1)/conv_D_enm #in Debye units

  ave,std=numpy.average(moments_normD),numpy.std(moments_normD)

  hist_spec_dipole,bin_edges=numpy.histogram(moments_normD,bins=1000,range=[0.0,10.0],density=True)
  bin_edges=bin_edges[:-1]
  bin_edges,hist_spec_dipole= bin_edges.reshape(-1,1),hist_spec_dipole.reshape(-1,1)
  hist_display = numpy.hstack((bin_edges,hist_spec_dipole))
  numpy.savetxt(outfile,hist_display)

'''
  for species in cdict:
    outfile=outstr+'_'+species+'.dat'
    spec_dipole = numpy.empty(0) #collection of all 1monomer D values.
    milist=midict[species] #residue index array
    for mindex in milist: #individual molecule.
      alist=topology.select('resid '+str(mindex))
      traj_1mon=traj.atom_slice(alist)
      moments=md.dipole_moments(traj_1mon,carray[alist])
      moments_normD=numpy.linalg.norm(moments,axis=1)/conv_D_enm #in Debye units
      spec_dipole=numpy.append(spec_dipole,moments_normD)
    print('Species {} : collected molecular dipole data {}'.format(species,len(spec_dipole)))
    ave,std=numpy.average(spec_dipole),numpy.std(spec_dipole)
    print('Species {} : Dipole ave {:8.3f} D std {:8.3f} D'.format(species,ave,std))

    hist_spec_dipole,bin_edges=numpy.histogram(spec_dipole,bins=1000,range=[0.0,10.0],density=True)
    bin_edges=bin_edges[:-1]
    bin_edges,hist_spec_dipole= bin_edges.reshape(-1,1),hist_spec_dipole.reshape(-1,1)
    hist_display = numpy.hstack((bin_edges,hist_spec_dipole))
    numpy.savetxt(outfile,hist_display)
'''
if __name__ == "__main__": main()

