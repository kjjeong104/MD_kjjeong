#!/home/kjeong23/softwares/bin/python3.4

import math
import sys
import numpy
import timeit
import mdtraj as md
from scipy.cluster import hierarchy
from scipy.spatial import distance 

def dmatrix3d(crd,box): #distance matrix in 3d
  for k in range(3):
    dmat1d=distance.pdist(crd[:,k].reshape(crd.shape[0],1))
    dmat1d[dmat1d>box[k]*0.5]-=box[k]
    try:
      dmat+=dmat1d**2
    except NameError:
      dmat=dmat1d**2
  #dmat[dmat>box*0.5]-=box
  dmat=numpy.sqrt(dmat)
  return dmat

def main():
  grofile=sys.argv[1]
  trjfile=sys.argv[2]
  dthres=0.85
  nthres=6

  start_time=timeit.default_timer()

  traj=md.load(trjfile,top=grofile)
  tailstr='C1'
  elapsed=timeit.default_timer() -start_time
  print('trajectory reading time passed {}'.format(elapsed))

  topology=traj.topology
  tends=topology.select('name '+tailstr)
  #for atom in topology.atoms:print(atom.name)
  #print(topology.atoms[0].name)
  nsurf=len(tends)
  nstep=traj.n_frames
  box=numpy.array(traj.unitcell_lengths[0])
  crd=traj.xyz
  totmic,onemic,already=[],[],[]

  for sindex in range(nstep):
    onestep_crd=crd[sindex]
    #currt=traj[sindex].time
    #print(currt)
    #print(tends,len(tends))
    te_crd=onestep_crd[tends]
    dmat=dmatrix3d(te_crd,box)
    cluslist=hierarchy.linkage(dmat,method='single') #single-linkage clustering
    efcllist=numpy.empty((0,4)) #effective cluster list
    #print(cluslist)
    for line in cluslist:
      if line[2]<=dthres: 
        #print(line)#distance criterion
        efcllist=numpy.vstack((efcllist,line)) 
      #if line[2]<=dthres: print(line)
    #for line in efcllist:
    #  print([int(line[0]),int(line[1]),line[2],int(line[3])])
    efcllist=efcllist[::-1] #flip matrix
    nlrow=efcllist.shape[0]
    #print(efcllist)
    for line in efcllist: #dividing cluster info into micelles
      onemic=[]
      templist=[int(line[0]),int(line[1])]
      while True:
        for x in templist:
          if x not in already:
            already.append(x)
            if x<nsurf: #no further linkage, sole molecule: include
              onemic.append(x)
              templist.remove(x)
            else:  #a minor cluster
              #old_line_ind=memnum-nsurf flip:nlrow-oldind-1
              flipind=nlrow-x+nsurf-1
              templist.append(int(efcllist[flipind][0]))
              templist.append(int(efcllist[flipind][1]))
          else:templist.remove(x)
        if len(templist)==0:break
      if len(onemic)>=nthres:
        onemic.sort()
        totmic.append(onemic)
    #print(dmat,numpy.size(dmat))
    #for line in totmic:print(line)
    print(len(totmic))
    print('aggN list')
    for row in totmic:
      print (len(row))

  elapsed=timeit.default_timer() -start_time
  print('clustering done {}'.format(elapsed))

  #print(nsurf,nstep,box)
  #print(crd[0][1][1])
  #print(traj.xyz)
  #print(tends)
  #print(box[0],box[1],box[2])

if __name__ == "__main__": main()
