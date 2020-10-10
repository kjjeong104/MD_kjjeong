#!/usr/bin/env python
# an intermediate output printer program from spatd_hbond_counter_v02 code
# creates hydrogen bond 'matrix' and print as an ASCII file.
# matrix file format: rows: snapshots. column: by molecule indices ( record the number of hydrogen bonds at that snapshot )

# program to calculate distance-angle distribution function of hbond trio.
# algorithm : read dcd file to get coord. also, get donor-hydrogen-acceptor atomtype series.
# ex) HBD molecule: urea (NH), HBA: Cl. then Nm Hm Cl
# determine 'upper bound distance' for hydrogen-acceptor distance (ex 5.0A) to check.
# also determine bin widths for r, cos(theta).
# collect D-A distances within threshold, then make angle calculation atom index list.
# calc cosine of dha angle. collect into 2-dimensional bins
#inputs: pdbfile, traj dcd file. interactive input for settings parameters
#not supporting 'spherical or ecliptical hydrogen bond cutoff definition' yet.
#also, use angstrom for basic length unit (be careful that mdtraj use nanometers as basic distance unit)
#v011 : use 'iterload' for efficient memory management.
#v011 : functional enhancement : not only counts donation, but also counts acception.

import math
import sys
import numpy
import timeit
import mdtraj as md
import copy

angtonm=0.10
degtorad=math.pi/180.0
#tstep=1.000 #used ps time unit. 1 snapshot : 1ps.

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  hbmatfile = sys.argv[3]
  #nhboutfile = sys.argv[4] #short time n_HB ave outfile
  #nsoutfile = sys.argv[5] #number of solid molecules outfile (time evolution)

  aname1=input("atom names for donors? ex) OW \n")
  aname2=input("atom names for hydrogens? ex) HW1 HW2 \n")
  aname3=input("atom names for acceptors? ex) UO \n")
  hbrcut=float(input("cutoff distance(D-A) for hydrogen bond in angstroms? ex) 3.5 ? \n"))*angtonm
#In the donor-H -- acceptor triplet, D - A RDF r_min should be considered.
  hbacut=float(input("cutoff angle in degree width from 180, for hydrogen bond? ex) 30 \n"))
  hbamin,hbamax=(180.0-hbacut)*degtorad,(180.0+hbacut)*degtorad #in radian.
  tskip=int(input("Once in how many frames do you want to take? ex) 10 \n"))
  teq=int(input("How many initial frames do you want to cut as equilibration? ex) 5000 \n"))
  #thbave=float(input("How many ps to calculate short-time n_HB average? ex) 20 \n"))
  #nhbscut=float(input("Desired short-time n_HB cutoff for solid identification? ex) 2.50 \n")) 
  nhbmax=4 #assume that a single molecule can H-bond 4 times at most
  start_time=timeit.default_timer()

  #loop of iterating chunkload of trajectory
  lchunk,flag=100,0
  trueistep,istep=0,0
  for fragtraj in md.iterload(trjfile,chunk=lchunk,top=topfile):
    if teq<=(trueistep+lchunk):
      fragtraj=fragtraj[::tskip]
      if flag==0: #initial setup
        topology=fragtraj.topology
        totnmon=topology.n_residues
        asplit1,asplit2,asplit3=aname1.split(),aname2.split(),aname3.split()
        text1,text2,text3='','',''
        for word in asplit1:
          text1+='name '+word+' or '
        for word in asplit2:
          text2+='name '+word+' or '
        for word in asplit3:
          text3+='name '+word+' or '
        text1,text2,text3=text1[:-4],text2[:-4],text3[:-4]
        seld,selh,sela=topology.select(text1),topology.select(text2),topology.select(text3)
        n_atomd,n_atomh,n_atoma=len(seld),len(selh),len(sela)
        print(n_atomd,n_atomh,n_atoma)
        dhpairs=[]
        for j in topology.bonds:
          if j[0].index in selh and j[1].index in seld:
            dhpairs.append([j[1].index,j[0].index])
          elif j[0].index in seld and j[1].index in selh:
            dhpairs.append([j[0].index,j[1].index])
        fulllist_angles,donmonindex_angles,accmonindex_angles=[],[],[]
        for row in dhpairs:
          for i in sela:
            if topology.atom(row[1]).residue!=topology.atom(i).residue:
              extrow=row.copy()
              extrow.append(i)
              fulllist_angles.append(extrow)
              extrow_moni=[topology.atom(x).residue.index for x in extrow]
              donmonindex_angles.append(extrow_moni[0])
              accmonindex_angles.append(extrow_moni[2])
        #list_dist=numpy.array(list_dist)
        fulllist_angles,donmonindex_angles=numpy.array(fulllist_angles),numpy.array(donmonindex_angles)
        uniqhbmon=numpy.unique(donmonindex_angles)
        nmon = len(uniqhbmon) #number of donor residues involved in the Hbond topology tracking.
        print(nmon,'hydrogen bond donating monomers ')
        n_angles_full = len(fulllist_angles)
        print("dhpairs # = {}, full list # angles = {} ".format(len(dhpairs),n_angles_full))
	#hbond counting information array:the whole
        allmon_hbond_count = numpy.zeros((totnmon,0),dtype=int)
        hbtot=0
        flag=1

      #passed the initial topology setup stage
      dist = (md.compute_distances(fragtraj,fulllist_angles[:,[0,2]])).flatten()
      angl = (md.compute_angles(fragtraj,fulllist_angles)).flatten()
      allmon_hbond_count1 = numpy.zeros((totnmon,lchunk),dtype=int)
      # recalculate NaN values of angles
      for nan_item in numpy.argwhere(numpy.isnan(angl)).reshape(-1):
        i_frame = int(nan_item/n_angles_full)
        i_angle = nan_item%n_angles_full
        #print(" Nan at {} th frame and {} th angle".format(i_frame,i_angle))
        #print(" <-- {} th atoms".format(list_angles[i_angle]))
        i_abc = fulllist_angles[i_angle]
        a = fragtraj.xyz[i_frame][i_abc[0]]
        b = fragtraj.xyz[i_frame][i_abc[1]]
        c = fragtraj.xyz[i_frame][i_abc[2]]
        #print(" <-- position: \n {} \n {} \n {}".format(a,b,c))
        ba = a - b
        bc = c - b
        cosine_angle = numpy.dot(ba, bc) / (numpy.linalg.norm(ba) * numpy.linalg.norm(bc))
        #print(" distance= {}".format(numpy.linalg.norm(bc)))
        angle = numpy.arccos(cosine_angle)
        #print(" get correct value from NaN, {} (rad) {} (deg)".format(angle, angle*180.0/numpy.pi))
        angl[nan_item] = copy.copy(angle)
      i_thresd,i_thresa=numpy.where(dist<=hbrcut),numpy.where(numpy.logical_and(hbamin<=angl, angl<=hbamax))
      i_hb = numpy.intersect1d(i_thresd[0],i_thresa[0])
      hbtot+=len(i_hb)
      #hbond monomer counting section.
      for x in i_hb:
        i_frame = int(x/n_angles_full) #local frame number. 
        i_angle = x%n_angles_full
        i_donmon,i_accmon = donmonindex_angles[i_angle],accmonindex_angles[i_angle]
        allmon_hbond_count1[i_donmon, i_frame] += 1
        allmon_hbond_count1[i_accmon, i_frame] += 1
      #end of 1 chunk. promote the indices, stack the data.
      trueistep+=lchunk
      istep+=lchunk
      allmon_hbond_count = numpy.hstack((allmon_hbond_count,allmon_hbond_count1))
      if istep%(10*lchunk)==0:
        elapsed=timeit.default_timer() - start_time
        print('finished snapshot {} time {}'.format(istep,elapsed))
     
    else: #skip.
      trueistep+=lchunk
  nstep=istep
  #average number of H-bonds normalization : when D-H -- A exist, number of A atoms satisfies criterion, per one D-H.
  #therefore, (Total count in the whole trajectory)/((# of frames)*(total# of D-H pairs in 1 system))
  hbavg=hbtot/(nstep*len(dhpairs))
  print('Detected H-bond rcut {:8.3f} A angcut {:8.3f} deg totcount {} avg {:11.4f}'.format(hbrcut/angtonm,hbacut,hbtot,hbavg))

  mon_hbond_count=allmon_hbond_count[uniqhbmon]
  hb_total_matrix=numpy.transpose(mon_hbond_count)
  numpy.savetxt(hbmatfile,hb_total_matrix,fmt='%d')

  elapsed=timeit.default_timer() - start_time
  print('finished job {}'.format(elapsed))

if __name__ == "__main__": main()

