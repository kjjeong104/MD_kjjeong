#!/usr/bin/env python
#program to read trajectory file, net charge value table and evaluate dipole moment 
#add conditional feature : use hydrogen bond criterion to sort out hydrogen bonds,
# and calculate conditional histogram , ave, std of dipole moment
# for molecule within the shell of Hbond. 
#inputs: pdbfile, traj (dcd or xtc) file, charge table file

#import math
import sys
import numpy
import mdtraj as md
import argparse
import timeit
import copy 

conv_D_enm = 0.020819434
angtonm=0.10
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
parser.add_argument('-o', '--outfile', default='dipole', nargs='?', 
  help='dipole moment histogram for the trajectory')
parser.add_argument('-b', '--begin',default=0,nargs='?',
  help='initial cut frames for equilibration ')
parser.add_argument('-skip', '--skip',default=1,nargs='?',
  help='frameskip ')
args = parser.parse_args()

def hbond_info_module(topfile):
  t=md.load(topfile)
  topology=t.topology
  aname1=input("atom names for donors? ex) OW \n")
  aname2=input("atom names for hydrogens? ex) HW1 HW2 \n")
  aname3=input("atom names for acceptors? ex) UO \n")
  hbrcut=float(input("cutoff distance(D-A) for hydrogen bond in angstroms? ex) 3.5 ? \n"))*angtonm
#In the donor-H -- acceptor triplet, D - A RDF r_min should be considered.
  hbacut=float(input("cutoff angle in degree width from 180, for hydrogen bond? ex) 30 or 180 \n"))
  hbamin,hbamax=(180.0-hbacut)*conv_deg_rad,(180.0+hbacut)*conv_deg_rad #in radian.

  #make atom indices list (before filtering too far pairs)
  #should avoid intramolecular atomic pair
  asplit1,asplit2,asplit3=aname1.split(),aname2.split(),aname3.split()
  text1,text2,text3='','',''
  for word in asplit1:
    text1+='name '+word+' or '
  for word in asplit2:
    text2+='name '+word+' or '
  for word in asplit3:
    text3+='name '+word+' or '
  text1,text2,text3=text1[:-4],text2[:-4],text3[:-4]
  seld=topology.select(text1)
  selh=topology.select(text2)
  sela=topology.select(text3)
  n_atomd,n_atomh,n_atoma=len(seld),len(selh),len(sela)
  print(n_atomd,n_atomh,n_atoma)
  dhpairs=[]
  for j in topology.bonds:
    if j[0].index in selh and j[1].index in seld:
      dhpairs.append([j[1].index,j[0].index])
    elif j[0].index in seld and j[1].index in selh:
      dhpairs.append([j[0].index,j[1].index])
  fulllist_angles,monindex_angles=[],[]
  for row in dhpairs:
    for i in sela:
      if topology.atom(row[1]).residue!=topology.atom(i).residue:
        extrow=row.copy()
        extrow.append(i)
        fulllist_angles.append(extrow)
        extrow_moni=[topology.atom(x).residue.index for x in extrow]
        monindex_angles.append(extrow_moni)
  #list_dist=numpy.array(list_dist)
  fulllist_angles,monindex_angles=numpy.array(fulllist_angles),numpy.array(monindex_angles)
  n_angles_full = len(fulllist_angles)
  print("dhpairs # = {}, full list # angles = {} ".format(len(dhpairs),n_angles_full))

  return fulllist_angles,monindex_angles,hbrcut,hbacut,hbamin,hbamax

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
    midict[rname]=numpy.append(midict[rname],int(i))
    carray=numpy.append(carray,cdict[rname])
  return carray,cdict,midict

#main fxn
def main():
  #Part to load coordinate file
  topfile = args.topfile
  trjfile = args.trjfile
  chgfile = args.chgfile
  outstr = args.outfile
  teq,tskip = int(args.begin),int(args.skip)
  #input 1 : load surf traj. (big file)
  fulllist_angles,monindex_angles,hbrcut,hbacut,hbamin,hbamax = hbond_info_module(topfile)
  n_angles_full=len(fulllist_angles)
  carray,cdict,midict=construct_carray_midict(topfile,chgfile)

  print('program started')
  start_time=timeit.default_timer()

  traj=md.load(trjfile,top=topfile)
  topology = traj.topology
  traj=traj[teq::tskip]
  nstep,nmon=traj.n_frames,topology.n_residues
  print(nstep, ' frames ',nmon,' molecules ')
  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))

  #prepare fragmental traj calc.(memory saver)
  if nstep>=200: #if there're many frames..
    nfrag=200
  else:
    nfrag=1
  #calculate individual molecular dipole moments.
  #For each species, collect differently.
  #now this part: should be combined with H-bond judgement algorithm.
  #step 1 : use H-bond judgement. For each monomer index, can make a 'digital matrix' to say
  #which MD frame is H-bonding, and which frame is not bonding. (this technique is useful in conditional rotcorr)
  mon_hbond_tf = numpy.zeros((nmon,nstep), dtype=bool) #boolean array of monomer hbond at each traj frame
  hbtot=0
  for ifrag in range(nfrag): #loop of fragmental calculation, to save memory.
    blength=int(nstep/nfrag)
    bstart,bend=ifrag*blength,(ifrag+1)*blength
    fragtraj=traj[bstart:bend]
    #calculate distances between donors and acceptors, angle 
    dist = (md.compute_distances(fragtraj,fulllist_angles[:,[0,2]])).flatten()
    angl = (md.compute_angles(fragtraj,fulllist_angles)).flatten()
    #print(" list # distances(within threshold) = {}".format(len(dist)))
    #print(" list # angles = {}".format(n_angles))
    # recalculate NaN values of angles
    for nan_item in numpy.argwhere(numpy.isnan(angl)).reshape(-1):
        i_frame = int(nan_item/n_angles_full)
        i_angle = nan_item%n_angles_full
        #print(" Nan at {} th frame and {} th angle".format(i_frame,i_angle))
        #print(" <-- {} th atoms".format(list_angles[i_angle]))
        i_abc = fulllist_angles[i_angle]
        a = traj.xyz[i_frame][i_abc[0]]
        b = traj.xyz[i_frame][i_abc[1]]
        c = traj.xyz[i_frame][i_abc[2]]
        print(" <-- position: \n {} \n {} \n {}".format(a,b,c))
        ba = a - b
        bc = c - b
        cosine_angle = numpy.dot(ba, bc) / (numpy.linalg.norm(ba) * numpy.linalg.norm(bc))
        print(" distance= {}".format(numpy.linalg.norm(bc)))
        angle = numpy.arccos(cosine_angle)
        print(" get correct value from NaN, {} (rad) {} (deg)".format(angle, angle*180.0/numpy.pi))
        angl[nan_item] = copy.copy(angle)
    i_thresd,i_thresa=numpy.where(dist<=hbrcut),numpy.where(numpy.logical_and(hbamin<=angl, angl<=hbamax))
    i_hb=numpy.intersect1d(i_thresd[0],i_thresa[0])
    hbtot+=len(i_hb)
    for x in i_hb: #hbond anglist array. For each frame, angle triplet indices are arranged.
      i_frame = int(x/n_angles_full) + bstart
      i_angle = x%n_angles_full
      i_mon = monindex_angles[i_angle]
      mon_hbond_tf[i_mon,i_frame] = True

    elapsed=timeit.default_timer() - start_time
    print('finished fragment {} H-bond calc. time {}'.format(ifrag,elapsed))
  print('Detected H-bond rcut {:8.3f} A angcut {:8.3f} deg totcount {}'.format(hbrcut/angtonm,hbacut,hbtot))

  #step 2 : get dipole moment data. apply boolean array to sort out non-Hbond frames for each monomer.
  #finalize dipole moment statistics.
  for species in cdict:
    outfile=outstr+'_'+species+'.dat'
    spec_dipole = numpy.empty(0) #collection of all 1monomer D values.
    milist=midict[species] #residue index array
    for mindex in milist: #individual molecule.
      mindex=int(mindex)
      alist=topology.select('resid '+str(mindex))
      traj_1mon=traj.atom_slice(alist)
      moments=md.dipole_moments(traj_1mon,carray[alist])
      moments_normD=numpy.linalg.norm(moments,axis=1)/conv_D_enm #in Debye units
      spec_dipole=numpy.append(spec_dipole,moments_normD[mon_hbond_tf[mindex]])
    print('Species {} : collected molecular dipole data {}'.format(species,len(spec_dipole)))
    ave,std=numpy.average(spec_dipole),numpy.std(spec_dipole)
    print('Species {} : Dipole ave {:8.3f} D std {:8.3f} D'.format(species,ave,std))

    hist_spec_dipole,bin_edges=numpy.histogram(spec_dipole,bins=1000,range=[0.0,10.0],density=True)
    bin_edges=bin_edges[:-1]
    bin_edges,hist_spec_dipole= bin_edges.reshape(-1,1),hist_spec_dipole.reshape(-1,1)
    hist_display = numpy.hstack((bin_edges,hist_spec_dipole))
    numpy.savetxt(outfile,hist_display)
    elapsed=timeit.default_timer() - start_time
    print('finished species {} dipole stat. time {}'.format(species,elapsed))

  #print(moments)
  #print(str(numpy.linalg.norm(moments[0])/conv_D_enm)+' D')

if __name__ == "__main__": main()

