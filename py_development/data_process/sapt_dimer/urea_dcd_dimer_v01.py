#!/home/kjeong23/softwares/bin/python3.4
# prototype program for urea dimer autoselection from MD traj.
# algorithm : read dcd file to get coord -> from given criterion, randomly pick dimers
# -> write xyz files

import math
import sys
import random
import numpy
import timeit
import mdtraj as md
from scipy.stats import itemfreq

#dimer search module. First pick a random monomer, then seek neighbor
def search_dimer(oneframe,crit):
  #pick random residue
  nmon,natoms,topology=oneframe.n_residues, oneframe.n_atoms, oneframe.topology
  m2=-1
  #if first guess fails, need to seek another. So make as while loop
  while True:
    m1=int(random.random()*nmon)
    atomids=topology.select('resid '+str(m1))
    search=[x for x in range(natoms) if x not in atomids]
    neigh=md.compute_neighbors(oneframe,crit[0],atomids,haystack_indices=search,periodic=False)
    resilist=[]
    for x in neigh[0]:
      resilist.append(topology.atom(x).residue.index)
    freq=itemfreq(resilist)
    for f in freq:
      if f[1]>=crit[1]:
        m2=f[0]
        break
    if m2!=-1:break  

  #print(m1,m2)
  dimer=oneframe.atom_slice(topology.select('resid '+str(m1)+' or resid '+str(m2)) )
  atomids=topology.select('resid '+str(m1))
  #dimer=dimer.image_molecules(inplace=False,anchor_molecules=atomids)
  return dimer

def main():
  #declare topology (.pdb) file and trajectory (.dcd) file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]

  #settings for dimer choice.
  ndim=int(input("how many dimer output files?\n"))
  cflag=int(input("criterion type? 1: atomic contact distance, 2: specific pair(unable), 3: H-bond(unable)\n"))
  if cflag==1:
    crit1=float(input("The threshold contact distance? (in A)\n"))
    crit1*=0.1 #unit conversion from angstrom to nanometers
    crit2=int(input("Number of required atom-atom contact for dimer selection?\n"))
    crit=[crit1,crit2]

  start_time=timeit.default_timer()
  #load files and prepare parameters
  oritopol=md.load(topfile).topology
  table,bonds=oritopol.to_dataframe()
  table.loc[ table['name'] == 'O2a','name' ] = 'O'
  table.loc[ table['name'] == 'Nm1','name' ], table.loc[ table['name'] == 'Nm2','name' ] = 'N','N'
  table.loc[ table['name'] == 'C2a','name' ] = 'C'
  table.loc[ table['name'] == 'Hm1','name' ], table.loc[ table['name'] == 'Hm2','name' ],\
  table.loc[ table['name'] == 'Hm3','name' ], table.loc[ table['name'] == 'Hm4','name' ] = 'H','H','H','H'
  topology=md.Topology.from_dataframe(table,bonds)
  traj=md.load_dcd(trjfile,top=topology)
  nstep=traj.n_frames
  traj=traj.atom_slice(topology.select('name != DO and name != DN1 and name != DN2 and name != DC'))
  
  #loop
  for dindex in range(ndim):
    #load 1 frame, then select dimer
    ntarf=int(nstep/2 + (nstep/2.0)/ndim * dindex) #index of target frame
    oneframe=traj[ntarf]
    #print(oneframe)
    sel_dimer=search_dimer(oneframe,crit)

    #write xyz file
    sel_dimer.save_xyz('urea_urea_'+str(dindex)+'.xyz')

    if dindex%50==0:
      elapsed=timeit.default_timer() - start_time
      print('step {} time {:10.4f}'.format(dindex,elapsed))

if __name__ == "__main__": main()
