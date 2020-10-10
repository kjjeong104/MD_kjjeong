#!/home/kjeong23/softwares/bin/python3.4
# prototype program for urea dimer autoselection from MD traj.
# algorithm : read dcd file to get coord -> from given criterion, randomly pick dimers
# -> write xyz files

import math
import sys
import numpy
import timeit
import mdtraj as md
from scipy.stats import itemfreq

def main():
  #declare topology (.pdb) file and trajectory (.dcd) file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  drudes=['DO','DN1','DN2','DC']
  #settings for dimer choice.
  #ndim=int(input("how many dimer output files?\n"))

  start_time=timeit.default_timer()

  #load files and process
  oritopol=md.load(topfile).topology
  table,bonds=oritopol.to_dataframe()
  table.loc[table['name'] == 'O2a','name'] = 'O'
  topology=md.Topology.from_dataframe(table,bonds)
  traj=md.load_dcd(trjfile,top=topology)
  
  #print(traj)
  traj=traj.atom_slice(topology.select('name != DO and name != DN1 and name != DN2 and name != DC'))
  oneframe=traj[100]
  topology=oneframe.topology
  atoms=topology.atoms
  #print(table.head())
  #print(oneframe)
  #print(topology)
  sel_dimer=oneframe.atom_slice(topology.select('resid 0 or resid 2'))
  #print(topology.select('resid 0 or resid 2'))
  #print(sel_dimer)
  #print(sel_dimer)
  #sel_dimer.image_molecules()
  sel_dimer.save_xyz('test.xyz')
  #apair=[[0,1]]
  #cont_result=md.compute_contacts(oneframe,contacts=apair,scheme='closest')
  #print(cont_result)
 
  natoms=oneframe.n_atoms
  atomids=topology.select('resid 0')
  print(atomids)
  search=[x for x in range(natoms) if x not in atomids]
  neigh=md.compute_neighbors(oneframe,0.3,atomids,haystack_indices=search,periodic=True)
  print(neigh,len(neigh[0]))

  resilist=[]
  for x in neigh[0]:
    resilist.append(topology.atom(x).residue.index)
  print(resilist)
  freq=itemfreq(resilist)
  #print(freq)
  for f in freq:
    if f[1]>=3:
      m2=f[0]
      break
  print(m2)
  #print(resilist)
  #loop
  #for dindex in range(ndim):

  #  if dindex%50==0:
  #    elapsed=timeit.default_timer() - start_time
  #    print(''.format(dindex,elapsed))
  elapsed=timeit.default_timer() - start_time
  print('elapsed time {}'.format(elapsed))

if __name__ == "__main__": main()
