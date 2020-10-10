#!/home/kjeong23/softwares/bin/python3.4
# prototype program for micelle detection from coordinate file
# algorithm: get coord&info -> starts loop. -> distance test(consider pbc)
# -> include -> further test again&again -> full 1micelle list -> store
# -> start looking for another micelle
# v05 : able to read trajectory ascii file for statistical treatment
# ** requires pbc -whole treatment
# v052 : gyration tensor eigenvalues display
# v061: binary data compatibility, v062 parallelization,
# also reconsideration on micelle judgement threshold decision

import math
import sys
import numpy
import timeit
import mdtraj as md
from scipy.cluster import hierarchy
from scipy.spatial import distance 

#dictionary for atomic/united atomic masses. V05.
amass={'C1':15.035, 'C2':14.027, 'C3': 14.027, 'C4': 14.027, 'C5': 14.027, \
'C6': 14.027, 'C7': 14.027, 'C8': 14.027, 'C9': 14.027, 'CA': 14.027, \
'P1': 30.9738, 'OD1': 15.9994, 'OD2': 15.9994, 'OD3': 15.9994}

def dmatrix3d(crd,box): #distance matrix in 3d
  for k in range(3):
    dmat1d=distance.pdist(crd[:,k].reshape(crd.shape[0],1))
    dmat1d[dmat1d>box[k]*0.5]-=box[k]
    try:
      dmat+=dmat1d**2
    except NameError:
      dmat=dmat1d**2
  dmat=numpy.sqrt(dmat)
  return dmat

def pbc_refine(tempcrd,box,acode): #coord refinement for wall-crossing micelle atoms
  move=[]
  natoms=tempcrd.shape[0]
  nmol=0
  tends=[]
  for i in range(natoms):
    if acode[i]=='C1':
      tends.append(i) #indices list for tail-end atoms
      nmol+=1
  namol=tends[1]-tends[0]
  for i in range(nmol):
    for j in range(i,nmol):
      for k in range(3):
        if (tempcrd[tends[j]][k]-tempcrd[tends[i]][k]) > (box[k]/2.0): #if j is cut by (1 boxlength) edge
          if [tends[j],k] not in move: move.append([tends[j],k])
        elif (tempcrd[tends[i]][k]-tempcrd[tends[j]][k]) > (box[k]/2.0): #if i is cut by (1 boxlength) edge
          if [tends[i],k] not in move: move.append([tends[i],k])

  for x in move:
    for y in range(namol):
      tempcrd[x[0]+y][x[1]]-=box[x[1]]    
  return tempcrd

def morphol(tempcrd,box,acode): #acode: list of atoms' names for a micelle
  #calculate morphology(COM,rg,asphericity) of an atom group. pbc refinement required a priori.
  #masses should be also considered by using dictionary!
##part 1 : calculating COM, also oc ###v051 performance upgrade attempt
  c,oc=numpy.zeros(3),numpy.zeros(3)
  natoms=tempcrd.shape[0]
  totmass,noxy=0.0,0
  oxys=[]
  for i in range(natoms):
    an=acode[i]
    m1=amass[an] #corresponding atom's mass
    c+=m1*tempcrd[i]
    totmass+=m1
    if an=='OD1' or an=='OD2' or an=='OD3':
      oc+=tempcrd[i]
      noxy+=1
      oxys.append(i)
  c/=totmass
  oc/=noxy
  #COM adjustment(put inside box)
  outbox=1
  while outbox==1:
    outbox=0
    for k in range(3):
      if c[k] < 0.0:
        c[k]+=box[k]
        outbox=1
      elif c[k] > box[k]:
        c[k]-=box[k]
        outbox=1
##part 2: micelle radius of gyration.(using only oxygens.)
##part 3: asphericity (using only oxygens. So, calculate oxygens' COM.)
 #definition T_mn = 1/2N^2 sum(i)sum(j)(r(i)_m-r(j)_m)(r(i)_n-r(j)_n)
 #or T_mn = 1/N sum(i) r(i)_m * r(i)_n when origin is set to COM
  Tg=numpy.zeros((3,3)) #gyration tensor. need to be 3*3 array
  for m in range(3): # dimension index 1
    for n in range(3): # dimension index 2
      for x in oxys:
        Tg[m,n]+=(tempcrd[x][m]-oc[m])*(tempcrd[x][n]-oc[n])
      Tg[m,n]/=noxy
  #diagonalize T (Also, Rg^2 = lmda_x^2 + lmda_y^2 + lmda_z^2)
  e_values, e_vectors = numpy.linalg.eig(Tg)
  e_values=numpy.sort(e_values) #sorting in ascending order
  rg2=e_values[0]+e_values[1]+e_values[2]
  rg=numpy.sqrt(rg2)
  e10,e20,e21 = e_values[1]-e_values[0],e_values[2]-e_values[0],e_values[2]-e_values[1]
  asp=(e10*e10+e20*e20+e21*e21)/(2.0*rg2*rg2) #asphericity

  Ti=numpy.zeros((3,3)) #moment of inertia tensor divided by N. get values from gyr tensor
  Ti[0,0],Ti[1,1],Ti[2,2] = Tg[1,1]+Tg[2,2], Tg[0,0]+Tg[2,2], Tg[0,0]+Tg[1,1] #diagonal elements
  Ti[0,1],Ti[0,2],Ti[1,0],Ti[1,2],Ti[2,0],Ti[2,1]=-Tg[0,1],-Tg[0,2],-Tg[1,0],-Tg[1,2],-Tg[2,0],-Tg[2,1] #off-diagonals
  ei_values, ei_vectors = numpy.linalg.eig(Ti) #diagonalize
  ei_values=numpy.sort(ei_values) #sorting in ascending order
  #semiaxis square lengths
  sv2=numpy.array([2.5*(ei_values[1]+ei_values[2]-ei_values[0]),\
                2.5*(ei_values[0]+ei_values[2]-ei_values[1]),\
                2.5*(ei_values[0]+ei_values[1]-ei_values[2])]) #sv2[0] is the largest
  #sv=numpy.sqrt(sv) #semiaxis vector (ref: J. Chem. Phys., Vol. 101, No. 10, 15 November 1994)
  
 #printing/returning section
  rg=round(rg,4)
  c=numpy.array([round(c[0],3),round(c[1],3),round(c[2],3)])
  sv2=numpy.array([round(sv2[0],4),round(sv2[1],4),round(sv2[2],4)])
  d=numpy.array([c[0],c[1],c[2],rg,asp,sv2[0],sv2[1],sv2[2]])
  return d

#main fxn:main administrator of algorithm progression
def main():
  #Part to load coordinate file. can process binary type. V06.
  grofile = sys.argv[1]
  trjfile = sys.argv[2] #now should load xtc file itself, not ASCII file
  outfile = open(sys.argv[3],'w') #output file for result summary
  comfile = open(sys.argv[4],'w') #output file for reduced expression of micelle
  dthres=float(input("threshold for micelle detection by tail position? recomm:0.85\n"))
  nthres=int(input("minimum aggregation number to be recognized as a micelle?recomm:6\n"))

  start_time=timeit.default_timer()
  surfstr,tailstr,namol='DEP','C1',14
  traj=md.load(trjfile,top=grofile)

  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))
  topology=traj.topology
  tends=topology.select('name '+tailstr)
  totacode=[]
  for atom in topology.atoms:totacode.append(atom.name)
  nsurf=len(tends) #length of C1 atoms' list
  nstep=traj.n_frames
  box=numpy.array(traj.unitcell_lengths[0])
  crd=traj.xyz
  totmic,onemic,already=[],[],[]

  #clustering section
  for sindex in range(nstep):
    onestep_crd=crd[sindex]
    te_crd=onestep_crd[tends]
    dmat=dmatrix3d(te_crd,box)
    cluslist=hierarchy.linkage(dmat,method='single') #single-linkage clustering
    efcllist=numpy.empty((0,4)) #effective cluster list
    for line in cluslist:
      if line[2]<=dthres: 
        efcllist=numpy.vstack((efcllist,line)) 
    efcllist=efcllist[::-1] #flip matrix
    nlrow=efcllist.shape[0]
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

    #com reduction, morphology calculation 
    coms,rgs,asps,sv2a,sv2b,sv2c=numpy.empty((0,3),float),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([])
    for row in totmic: #tailend member list of 1 micelle
      mic_atoms_list=numpy.array([])
      acode=[]
      for i in row:
        temp_list=topology.select('resid '+str(i))
        mic_atoms_list=numpy.concatenate([mic_atoms_list,temp_list])
      mic_atoms_crd=numpy.empty((0,3))
      for i in mic_atoms_list:
        mic_atoms_crd=numpy.vstack( (mic_atoms_crd,onestep_crd[i]) )
        acode.append(totacode[int(i)])
      refined_mic_crd=pbc_refine(mic_atoms_crd,box,acode)
      elapsed=timeit.default_timer() - start_time
      #print('debug:micelle pbc refinement {}'.format(elapsed))
      d=morphol(refined_mic_crd,box,acode) #instead of total data array, put 1 micelle atoms crd array
      #c=md.compute_center_of_mass(mic_atoms_crd)
      c=numpy.array([d[0],d[1],d[2]])
      coms=numpy.vstack((coms,c))
      #rg=md.compute_rg(mic_atoms_crd)
      rgs=numpy.append(rgs,d[3])
      asps=numpy.append(asps,d[4])
      sv2a,sv2b,sv2c=numpy.append(sv2a,d[5]),numpy.append(sv2b,d[6]),numpy.append(sv2c,d[7]) 

    #printing section
    micindex=0
    outfile.write('step# {} numofmcls {}\n'.format(sindex,len(totmic)))
    for row in totmic:
      micindex+=1
      outfile.write('mic# {:3} N= {:3} r= {:8.4f} A= {:8.4f} sv2abc= {:8.4f} {:8.4f} {:8.4f} members: {}\n'.\
format(micindex,len(row),rgs[micindex-1],round(asps[micindex-1],4),sv2a[micindex-1],sv2b[micindex-1],sv2c[micindex-1],row))
    #outfile.write('{} Unaggregated molecules: {}\n'.format(unagg,hom))

    #com output writing section
    outindex=0
    comfile.write('Centers of micelles: step {}\n'.format(sindex))
    comfile.write(' {}\n'.format(micindex))
    for row in coms:
      outindex+=1
      comfile.write('{:5}{:5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}\n'.format(outindex,'COM','X',outindex,row[0],row[1],row[2]))
    comfile.write('      {:8.3f}  {:8.3f}  {:8.3f}\n'.format(box[0],box[1],box[2]))

    elapsed=timeit.default_timer() - start_time
    print('finished calculating step {} time {}'.format(sindex,elapsed))
 
  outfile.close()
  comfile.close()

if __name__ == "__main__": main()

