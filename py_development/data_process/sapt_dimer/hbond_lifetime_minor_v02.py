#!/home/kjeong23/softwares/bin/python3.4
# program to calculate hydrogen bonding lifetime.
# algorithm : read dcd file to get coord. also, get donor-hydrogen-acceptor atomtype series.
# for specific d-h-a triplet, get pre-determined criterion for h-a distance, d-h-a angle.
# ex) HBD molecule: urea (NH), HBA: Cl. then Nm Hm Cl
# at first, analyze topology to get list of all available hydrogen bonds of designated type.
# for every snapshot, analyze all combinations(necessary!) to get list of hydrogen bonds.
# should collect <h> values from each snapshot, to calculate <h> of the entire ensemble.
# (h=1 for H-bond, h=0 for non-Hbond for a specific atom triplet)
# By comparing the list of hydrogen bonds of each snapshot, calculate h(t)
# C(t) = (<h(0)h(t)> - <h>^2) / (<h2>-<h>2)
# citation : J. Chem. Phys. 148, 193843 (2018)
#inputs: pdbfile, traj dcd file
#output: file of (t, C(t) columns) for hydrogen bond time correlation function
# be careful that mdtraj use nanometers as basic distance unit
# be able to separate topologies in different cpus
# minor version: after getting computation time issue, just try to first print
# H-bond topology index for each snapshot.

import math
import sys
import numpy
import timeit
import mdtraj
import copy

#tstep=1.000 #used ps time unit. 1 snapshot : 1ps.

#def autocorr_simple(x): # <x(0)x(t)> without considering <x>,<x2>
#  tlimit=int(len(x)/2)
  #a,s2=numpy.mean(x),numpy.var(x)
#  zerot,zerotcount=numpy.zeros(tlimit),numpy.zeros(tlimit)
  #zerot[0],zerotcount[0]=1.000,1.000
#  for j in range(0,tlimit):
#    for i in range(len(x)-j):
#      zerot[j]+=x[i]*x[i+j]
#      zerotcount[j]+=1
  #ct=numpy.divide(ct,ctcount)

#  return zerot,zerotcount

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  outfile = open(sys.argv[3],'w')

  aname1=input("atom names for donors? ex) Nm1 Nm2 \n")
  aname2=input("atom names for hydrogens? ex) Hm1 Hm2 Hm3 Hm4 \n")
  aname3=input("atom names for acceptors? ex) CL \n")
  rrange=input("minimum, maximum H-A distance in nanometers, to judge as hydrogen bond? ex) 0.00 0.30 \n")
  arange=input("angle range in cosines, to judge as hydrogen bond? ex) -1.0 -0.6\n ")
  rsplit,ansplit=rrange.split(),arange.split()
  rmin,rmax,amin,amax=float(rsplit[0]),float(rsplit[1]),float(ansplit[0]),float(ansplit[1])
  #mi12=input("What is the residue number index interval of MSD calculation? ex) 0 299 \n")
  #mi12=mi12.split()
  #mi1,mi2=int(mi12[0]),int(mi12[1])
  tskip=int(input("Once in how many frames do you want to take? ex) 1 \n"))
  teq=int(input("How many initial frames do you want to cut as equilibration? ex) 5000 \n"))
  #tstep=float(input("How much is the timstep between snapshots in dcd file (in ps)? ex) 10\n"))
  #tstep=tstep*float(tskip)
  #toprange=input("D-H--A triplet index range for multi-cpu calc? ex) 0 10000 (if unnecessary put -1 -1) \n")
  #if toprange=="-1 -1":
  #  toprange1,toprange2=-1,-1
  #else:
  #  toprsplit=toprange.split()
  #  toprange1,toprange2=int(toprsplit[0]),int(toprsplit[1])

  start_time=timeit.default_timer()

  #input 1 : load surf traj. (big file)
  traj=mdtraj.load(trjfile,top=topfile)
  traj=traj[teq::tskip]
  topology=traj.topology
  nstep=traj.n_frames
  nmon=topology.n_residues
  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))
  print(nstep,nmon)

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
  fulllist_angles=[]
  for row in dhpairs:
    for i in sela:
      if topology.atom(row[1]).residue!=topology.atom(i).residue:
        extrow=row.copy()
        extrow.append(i)
        fulllist_angles.append(extrow)
  #list_dist=numpy.array(list_dist)
  fulllist_angles=numpy.array(fulllist_angles)
  #if toprange1!=-1:
  #  if toprange2>len(fulllist_angles):
  #    toprange2=len(fulllist_angles)
  #  fulllist_angles=fulllist_angles[toprange1:toprange2]
  n_angles_full = len(fulllist_angles)
  print(" list # angles = {} ".format(n_angles_full))

  #should construct h-matrix (h values for 1 snapshot * all snapshots)
  #but the full h-matrix is too large, consuming too much memory. (ex) 320000*2000 int values)
  #so, slice the topology into smaller dha triplets. Different dha triplets' h-values are uncorrelated.
  #if n_angles_full > 10000:
  #  topbl=10000
  #  ntopb=int(numpy.ceil(float(n_angles_full)/topbl))
  #else:
  #  topbl,ntopb=n_angles_full,1

  #Do not slice topology. Instead, just print H-bonded topology number.
  for i_frame in range(nstep):
    traj1=traj[i_frame]
    dist1=(mdtraj.compute_distances(traj1,fulllist_angles[:,[1,2]])).flatten()
    angl1=(mdtraj.compute_angles(traj1,fulllist_angles)).flatten()
    for nan_item in numpy.argwhere(numpy.isnan(angl1)).reshape(-1):
      i_angle = nan_item%len(angl1)
      i_abc = fulllist_angles[i_angle]
      print(i_abc)
      print(traj1.xyz.shape)
      #print(traj1.xyz)
      a = traj1.xyz[0][i_abc[0]]
      b = traj1.xyz[0][i_abc[1]]
      c = traj1.xyz[0][i_abc[2]]
      print(" <-- position: \n {} \n {} \n {}".format(a,b,c))
      ba = a - b
      bc = c - b
      cosine_angle = numpy.dot(ba, bc) / (numpy.linalg.norm(ba) * numpy.linalg.norm(bc))
      print(" distance= {}".format(numpy.linalg.norm(bc)))
      angle = numpy.arccos(cosine_angle)
      print(" get correct value from NaN, {} (rad) {} (deg)".format(angle, angle*180.0/numpy.pi))
      angl1[nan_item] = copy.copy(angle)
    cosangl1=numpy.cos(angl1)
    index_thresd=numpy.where(numpy.logical_and(rmin<=dist1, dist1<=rmax))
    index_thresa=numpy.where(numpy.logical_and(amin<=cosangl1, cosangl1<=amax))
    index_thresda=numpy.intersect1d(index_thresd,index_thresa)

    #direct output printing
    #numpy.set_printoptions(linewidth=numpy.inf)
    outstr=numpy.array2string(index_thresda, max_line_width=10000)
    outstr=outstr.replace(",","")
    outstr=outstr.replace("[","")
    outstr=outstr.replace("]","")
    outfile.write("{}\n".format(outstr))

    if i_frame%20==0:
      elapsed=timeit.default_timer() - start_time
      print('step {} time {}'.format(i_frame,elapsed))

  outfile.close()

if __name__ == "__main__": main()

