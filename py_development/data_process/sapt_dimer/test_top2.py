#!/home/kjeong23/softwares/bin/python3.4

import math
import sys
import numpy
import timeit
import mdtraj as md

#main fxn
def main():
  start_time=timeit.default_timer()
  #Part to load coordinate file
  topfile = sys.argv[1]
  #aname1=input("atom names for donors? ex) Nm1 Nm2 \n")
  #aname2=input("atom names for hydrogens? ex) Hm1 Hm2 Hm3 Hm4 \n")
  #aname3=input("atom names for acceptors? ex) CL \n")
  aname1,aname2,aname3='Nm1 Nm2','Hm1 Hm2 Hm3 Hm4','CL'
  #input 1 : load surf traj. (big file)
  topology=md.load(topfile).topology
  nmon=topology.n_residues

  #make atom indices list (before filtering too far pairs)
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
  #print(seld)
  #print(selh)
  #print(sela)
  print(n_atomd,n_atomh,n_atoma)
  dhpairs=[]

  for j in topology.bonds:
    if j[0].index in selh and j[1].index in seld:
      dhpairs.append([j[1].index,j[0].index])
    elif j[0].index in seld and j[1].index in selh:
      dhpairs.append([j[0].index,j[1].index])

  print(len(dhpairs))
  #print(dhpairs)
 # for i in topology.bonds:
 #   print(i[0].index)
  fulllist_angle=[]
  for row in dhpairs:
    for i in sela:
      extrow=row.copy()
      extrow.append(i)
      fulllist_angle.append(extrow)
  #list_dist=numpy.array(list_dist)
  #print(fulllist_angle)
  fulllist_angle=numpy.array(fulllist_angle)
  n_angles_full = len(fulllist_angle)
  print(fulllist_angle)
  print(" full list # angles = {} ".format(n_angles_full))


  #list_dist=numpy.array(list_dist)
  #calculate distances between hydrogen and acceptors, filter them 
  # calculate distance

  #printing section - should regard gnuplot pm3d-compatible format.

  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed)) 

if __name__ == "__main__": main()

