#!/home/kjeong23/softwares/bin/python3.4
#Self-intermediate scattering function (F_s(q,t)) calculator from gromacs traj. Isotropic version.
# algorithm : read gro and traj. according to user input, 
#use the same loop structure and coordinate auto-unwrapping method as we calculate MSD, to calculate
# the displacement and following calc of imaginary exp.
# Assume that all atoms are equivalent, and we don't care about COM of molecules.
# input : grofile trjfile. output : Fsqt file.
# Formulas are listed below.
#direct calc: Fs(q,t)=1/N <sum(j=1~N) exp [-iq dot (r_j(t) - r_j(0))]> 
#FT of van Hove function : F(k,t)=integral G(r,t) exp(-ik dot r) dr
#Isotropic system : q dot r reduces to qr cos(theta), so the ens.ave reduces to
# sin(kr)/kr when r=r_j(t)-r_j(0). (r^2 factor doesn't exist in this case, because
# we only collect wrt orientations, and not swiping through r=-inf~inf.
# choice of q: regard qa=2pi when a is real space size of interest. Or can take 1st peak |q| in S(q)

import math
import sys
import numpy
import timeit
import mdtraj as md

def pbcdr(ri,rf,box): #displacement vector with considering pbc
  dr=rf-ri
  for i in range(3):
    if dr[i]>(box[i]/2.0):
      dr[i]-=box[i]
    elif dr[i]<(-box[i]/2.0):
      dr[i]+=box[i]
  return dr

def iso_selfint_calc(trimtraj,q,ntbin): #isotropic Fsqt calculator. assumes trajectory has 1 particle.
  fsqt1,count1=numpy.zeros(ntbin),numpy.zeros(ntbin)
  trimstep=trimtraj.n_frames
  crd=trimtraj.xyz
  box=numpy.array(trimtraj.unitcell_lengths[0])
  #coordinate unwrapping concern : automatically bypass unwrapping concern, by tracking 1 stepwise.
  for i in range(trimstep-1): #initial time
    vec=numpy.zeros(3)
    for j in range(i+1,i+1+ntbin):
      if i+1+ntbin<=trimstep:
        dr=pbcdr(crd[j][0],crd[j-1][0],box)
        vec+=dr
        rmag=numpy.linalg.norm(vec)
        qr=q*rmag
        if qr!=0: isoft=numpy.sin(qr)/(qr)
        else: isoft=1.0
        fsqt1[j-i-1]+=isoft
        count1[j-i-1]+=1
  return fsqt1,count1

def main():
  grofile = sys.argv[1]
  trjfile = sys.argv[2]
  outfile = open(sys.argv[3],'w')

  #user manual input info
  q=float(input("Value of |q| for Fs(q,t) calc (in nm-1)? ex) 1.80 \n"))
  tinfo=input("snapshot timestep(in ns), trange max (in #frames)? ex) 0.2 1000\n")
  tinfo=tinfo.split()
  tstep,tbin,tmin,tmax=float(tinfo[0]),1,1,int(tinfo[2])
  ntbin=(tmax-tmin)/tbin+1
  #ntbin=int(round(ntbin))
  mi12=input("What is the surfactant index interval of MSD calculation? ex) 0 9 \n")
  mi12=mi12.split()
  mi1,mi2=int(mi12[0]),int(mi12[1])

  start_time=timeit.default_timer()
 
  #input file loading
  traj=md.load(trjfile,top=grofile)
  topology=traj.topology
  if topology.n_residues <= mi2: #surfactant index exceeded
    mi2=topology.n_residues-1 #autofix surfactant index range
  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))
  print(traj)
  nstep=traj.n_frames

  #prepare bins
  fsqt,count=numpy.zeros(ntbin),numpy.zeros(ntbin)

  #calculation loop
  for mindex in range(mi1,mi2+1):
    trimtraj=trimtraj.atom_slice(topology.select('resid '+str(mindex)))
    fsqt1,count1=iso_selfint_calc(trimtraj,q,ntbin) #Fsqt calc for 1 atom (ens.ave)
    fsqt+=fsqt1
    count+=count1
    elapsed=timeit.default_timer() - start_time
    print('surf# {} Fs(q,t) calculated. time {:11.4f}'.format(mindex,elapsed))

  #printing section
  fsqt=numpy.divide(fsqt,count,out=numpy.zeros_like(fsqt), where=count!=0) #normalize
  #this code's bin doesn't include t=0. manually add Fs(q,t)=1 for t=0
  outfile.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(0,1,1))
  for i in range(ntbin):
    outfile.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(tstep*(i+tmin),fsqt[i],count[i]))
  outfile.close()

if __name__ == "__main__": main()
