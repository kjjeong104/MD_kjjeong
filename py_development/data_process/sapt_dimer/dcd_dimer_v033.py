#!/home/kjeong23/softwares/bin/python3.4
# program for dimer autoselection from MD traj.
# algorithm : read dcd file to get coord -> from given criterion, randomly pick dimers
# -> write xyz files
# v02 : generalize for arbitrary heterodimer. Also, considers vdw radius of atoms
# v021 : remove drude particles.
# v03 : added specific pair choice
# v032 : can remove other types of ghost sites. now in atomfile, set element name
# of drude sites as 'del' or 'delete'
# v033 : vdw contact surface introduced.
### Warning : This version allows internal atomic pair (For Fetisov et al. trajectory)
#Also, due to mixed residue order of Fetisov et al., trajectory, to avoid code crash,
#remove the part that substituting atom names

_ATOMIC_RADII = {'H'   : 0.120, 'He'  : 0.140, 'Li'  : 0.076, 'Be' : 0.059, 
                 'B'   : 0.192, 'C'   : 0.170, 'N'   : 0.155, 'O'  : 0.152, 
                 'F'   : 0.147, 'Ne'  : 0.154, 'Na'  : 0.102, 'Mg' : 0.086, 
                 'Al'  : 0.184, 'Si'  : 0.210, 'P'   : 0.180, 'S'  : 0.180, 
                 'Cl'  : 0.181, 'Ar'  : 0.188, 'K'   : 0.138, 'Ca' : 0.114, 
                 'Sc'  : 0.211, 'Ti'  : 0.200, 'V'   : 0.200, 'Cr' : 0.200, 
                 'Mn'  : 0.200, 'Fe'  : 0.200, 'Co'  : 0.200, 'Ni' : 0.163, 
                 'Cu'  : 0.140, 'Zn'  : 0.139, 'Ga'  : 0.187, 'Ge' : 0.211, 
                 'As'  : 0.185, 'Se'  : 0.190, 'Br'  : 0.185, 'Kr' : 0.202, 
                 'Rb'  : 0.303, 'Sr'  : 0.249, 'Y'   : 0.200, 'Zr' : 0.200, 
                 'Nb'  : 0.200, 'Mo'  : 0.200, 'Tc'  : 0.200, 'Ru' : 0.200, 
                 'Rh'  : 0.200, 'Pd'  : 0.163, 'Ag'  : 0.172, 'Cd' : 0.158, 
                 'In'  : 0.193, 'Sn'  : 0.217, 'Sb'  : 0.206, 'Te' : 0.206, 
                 'I'   : 0.198, 'Xe'  : 0.216, 'Cs'  : 0.167, 'Ba' : 0.149, 
                 'La'  : 0.200, 'Ce'  : 0.200, 'Pr'  : 0.200, 'Nd' : 0.200, 
                 'Pm'  : 0.200, 'Sm'  : 0.200, 'Eu'  : 0.200, 'Gd' : 0.200, 
                 'Tb'  : 0.200, 'Dy'  : 0.200, 'Ho'  : 0.200, 'Er' : 0.200, 
                 'Tm'  : 0.200, 'Yb'  : 0.200, 'Lu'  : 0.200, 'Hf' : 0.200, 
                 'Ta'  : 0.200, 'W'   : 0.200, 'Re'  : 0.200, 'Os' : 0.200, 
                 'Ir'  : 0.200, 'Pt'  : 0.175, 'Au'  : 0.166, 'Hg' : 0.155, 
                 'Tl'  : 0.196, 'Pb'  : 0.202, 'Bi'  : 0.207, 'Po' : 0.197, 
                 'At'  : 0.202, 'Rn'  : 0.220, 'Fr'  : 0.348, 'Ra' : 0.283, 
                 'Ac'  : 0.200, 'Th'  : 0.200, 'Pa'  : 0.200, 'U'  : 0.186, 
                 'Np'  : 0.200, 'Pu'  : 0.200, 'Am'  : 0.200, 'Cm' : 0.200, 
                 'Bk'  : 0.200, 'Cf'  : 0.200, 'Es'  : 0.200, 'Fm' : 0.200, 
                 'Md'  : 0.200, 'No'  : 0.200, 'Lr'  : 0.200, 'Rf' : 0.200, 
                 'Db'  : 0.200, 'Sg'  : 0.200, 'Bh'  : 0.200, 'Hs' : 0.200, 
                 'Mt'  : 0.200, 'Ds'  : 0.200, 'Rg'  : 0.200, 'Cn' : 0.200, 
                 'Uut' : 0.200, 'Fl'  : 0.200, 'Uup' : 0.200, 'Lv' : 0.200, 
                 'Uus' : 0.200, 'Uuo' : 0.200} 

import math
import sys
import random
import numpy
import timeit
import mdtraj as md
from scipy.stats import itemfreq

#dimer search module. First pick a random monomer, then seek neighbor
def search_heterodimer(oneframe,crit,rid1,rid2,cflag,ap1list,ap2list):
  #pick random residue
  topology=oneframe.topology
  nmon1=rid1[1]-rid1[0]+1
  nmon2=rid2[1]-rid2[0]+1
  m2=-1
  search=topology.select('resid '+str(rid2[0])+' to '+str(rid2[1])) #search target : monomer type 2 atoms
  if cflag==2: #if specific pair is required
    search=[val for val in search if val in ap2list]
  #if first guess fails, need to seek another. So make as while loop
  while True:
    m1=int(random.random()*nmon1)+rid1[0]
    atomids=topology.select('resid '+str(m1))
    if cflag==2:
      atomids=[val for val in atomids if val in ap1list]
    #search=[x for x in range(natoms) if x not in atomids] #search target : monomer type 2 atoms
    neigh=md.compute_neighbors(oneframe,crit[0],atomids,haystack_indices=search,periodic=True)
    resilist=[]
    for x in neigh[0]:
      resilist.append(topology.atom(x).residue.index)
    freq=itemfreq(resilist)
    for f in freq:
      if f[1]>=crit[1] and m1!=f[0]:
      #if f[1]>=crit[1]:
        m2=f[0]
        break
    if m2!=-1:break  

  #print(m1,m2)
  dimer=oneframe.atom_slice(topology.select('resid '+str(m1)+' or resid '+str(m2)) )
  atomids=topology.select('resid '+str(m1))
  #dimer=dimer.image_molecules(inplace=False,anchor_molecules=atomids)
  dimer=dimer.center_coordinates()
  return dimer

def main():
  #declare topology (.pdb) file and trajectory (.dcd) file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  #atomfile = open(sys.argv[3],'r')

  #settings for dimer choice.
  ndim=int(input("how many dimer output files?\n"))
  cflag=int(input("criterion type? 1: atomic contact distance, 2: specific pair, 3: vdw contact surface\n"))
  if cflag==1:
    crit1=float(input("The threshold contact distance? (in A)\n"))
    crit1*=0.1 #unit conversion from angstrom to nanometers
    crit2=int(input("Number of required atom-atom contact for dimer selection?\n"))
    crit=[crit1,crit2]
    ap1,ap2='None','None'
  elif cflag==2:
    apair=input("Atom names (in pdb file) for pairs? ex) H13 O2a\n")
    asplit=apair.split()
    ap1,ap2=asplit[0],asplit[1]
    crit1=float(input("The threshold contact distance? (in A)\n"))
    crit1*=0.1 #unit conversion from angstrom to nanometers
    crit2=int(input("Number of required atom-atom contact for dimer selection? ex) 1\n"))
    crit=[crit1,crit2]
  elif cflag==3:
    crit1=float(input("\n"))
    
    
  mon12=input("Name of two monomers? ex) EMIM Cl\n")
  msplit=mon12.split()
  namemon1,namemon2=msplit[0],msplit[1]
  rid_total=input("residue id bound for mon1 and mon2? ex) 0 149 150 299\n")
  rsplit=rid_total.split()
  rid1,rid2=numpy.array([int(rsplit[0]),int(rsplit[1])]),numpy.array([int(rsplit[2]),int(rsplit[3])])
  teq=int(input("How many initial frames do you want to cut as equilibration? ex) 5000 \n"))
  save_drude=input("Include drude particles? y or n \n")

  start_time=timeit.default_timer()

  #load atomfile (pdb name, actual element name, vdw radii)
  #pdbname1,pdbname2,element1,element2=[],[],[],[]
  #vdwr1,vdwr2=numpy.empty(0),numpy.empty(0)
  #flag=0
  #for aline in atomfile:
  #  if 'mon1' in aline:
  #    flag=1
  #  elif 'mon2' in aline:
  #    flag=2
  #  else:
  #    asplit=aline.split()
  #    if flag==1:
  #      pdbname1.append(asplit[0])
  #      element1.append(asplit[1])
  #    elif flag==2:
  #      pdbname2.append(asplit[0])
  #      element2.append(asplit[1])
  #natom1,natom2=len(pdbname1),len(pdbname2)
  #construct vdw hetero diameter matrix
  #not implemented yet

  #load files and prepare parameters
  #caution : can have bug if drude particle indices are preceding
  #oritopol=md.load(topfile).topology
  #oritopol=oritopol.subset([atom.index for atom in oritopol.atoms if ('D' not in atom.name)])
  #table,bonds=oritopol.to_dataframe()
  #for i in range(natom1):
  #  table.loc[ table['name'] == pdbname1[i],'name' ] = element1[i]
  #for i in range(natom2):
  #  table.loc[ table['name'] == pdbname2[i],'name' ] = element2[i]
  #topology=md.Topology.from_dataframe(table,bonds)
  topology=md.load(topfile).topology
  traj=md.load_dcd(trjfile,top=topology)
  traj=traj[teq:]
  nstep=traj.n_frames
  #drude particle exclusion
  if save_drude=='n':
    traj=traj.atom_slice([atom.index for atom in topology.atoms if ('del' not in atom.name) ])
    #ap1list,ap2list construction (reload topfile)
    #oritopol=oritopol.subset([atom.index for atom in oritopol.atoms if ('del' not in atom.name)])
  #ap1list,ap2list=oritopol.select('name '+ap1),oritopol.select('name '+ap2)
  ap1list,ap2list=topology.select('name '+ap1),topology.select('name '+ap2)
  elapsed=timeit.default_timer() - start_time
  print('traj loading complete. time {:10.4f}'.format(elapsed))
  print('nstep: {}'.format(nstep))
 
  #loop
  for dindex in range(ndim):
    #load 1 frame, then select dimer
    ntarf=int(nstep/ndim * dindex) #index of target frame
    oneframe=traj[ntarf]
    #print(oneframe)
    sel_dimer=search_heterodimer(oneframe,crit,rid1,rid2,cflag,ap1list,ap2list)

    #write xyz file
    sel_dimer.save_xyz(str(namemon1)+'_'+str(namemon2)+'_'+str(dindex)+'.xyz')
    sel_dimer.save_pdb(str(namemon1)+'_'+str(namemon2)+'_'+str(dindex)+'.pdb')

    if dindex%10==0:
      elapsed=timeit.default_timer() - start_time
      print('step {} time {:10.4f}'.format(dindex,elapsed))
 
  atomfile.close()

if __name__ == "__main__": main()
