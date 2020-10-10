#!/home/kjeong23/softwares/bin/python3.4
# program to calculate hoppingless MSD of surfactants
# v04 renovation: imported mdtraj to process xtc file. Also, xtc file is truncated with 1 atomtype.
# definition of non-hopping: until it is detected in 'other' micelle.
# in migevent scheme, until the arrival time for the destination micelle.
# algorithm: read migevent to record cutting points -> track 1 surfactant. cut traj.
# -> record r0 (now don't calc COM) -> calc square displacement until it hops.(consider PBC)
# in every step of COM SD calc, pile up in the bins, adding up "counts"(w.r.t. t)
# to compare statistical weight of all dr^2(t). -> calc MSD(t).
# -> display MSD(t), count(t).
# ** binary .xtc trajectory file: requires pbc -whole treatment!!
#inputs: grofile, traj xtc file, migfile (listed inter-micellar exchanges, molecule index sorted)
#output: file of {t, nonhop_MSD(t), count(t)}
#v03: started consideration of hopping: when a surfactant hopps, set new r0,t=0
#philosophy of ensemble averaging : the issue is, because of unequal length of non-hop traj,
#averaging without care will result in unequal weighting of averaging. but, still do that.
#leave statistical count beside, to make good averaging of different surfactant groups later.

import math
import sys
import numpy
import timeit
import mdtraj as md

tstep=0.200 #use skip10-trajectory

#now do not use atomic mass matrix anymore. Don't care about COM, and only track 1 atom.

#def mig_stepteller(migline,tb1,tb2,mi1,mi2): #read migline, record surf index and time
#  #read surf index and initial time(important)
#  split=migline.split()
#  mindex,tarr=int(split[0]),float(split[4]) #mindex : 'string' form
#  hop_info=[]
#  if mi1<=mindex<=mi2 and tb1<=tarr<=tb2:  #collect info only in this case
#    ntarr=(tarr-tb1)/tstep #collect info only in this case. frame index in entire trajectory
#    ntarr=int(round(ntarr))
#    hop_info=[mindex,ntarr]
#  return hop_info

def pbcdr(ri,rf,box): #displacement vector with considering pbc
  dr=rf-ri
  for i in range(3):
    if dr[i]>(box[i]/2.0):
      dr[i]-=box[i]
    elif dr[i]<(-box[i]/2.0):
      dr[i]+=box[i]
  return dr

def SD_calculator(trimtraj,nstep): #square displacement calculator. assumes that trajectory has 1 particle.
  sdbin1,count1=numpy.zeros(nstep),numpy.zeros(nstep)
  trimstep=trimtraj.n_frames
  crd=trimtraj.xyz
  box=numpy.array(trimtraj.unitcell_lengths[0])
  #for ntau in range(trimstep): #use tau-step loop first
  #  for i in range(trimstep-ntau): #initial snapshot
  for i in range(trimstep-1):
    vec=numpy.zeros(3)
    for j in range(i+1,trimstep):
      dr=pbcdr(crd[j][0],crd[j-1][0],box)
      vec+=dr
      rmag2=numpy.dot(vec,vec)
      sdbin1[j-i]+=rmag2
      count1[j-i]+=1
  sdbin1[0],count1[0]=0.0,1
  return sdbin1,count1

#main fxn
def main():
  #Part to load coordinate file, migration event file
  grofile = sys.argv[1]
  trjfile = sys.argv[2]
  migfile = sys.argv[3]
  outfile = open(sys.argv[4],'w')

  tblock=input("What is the time interval of trajectory you're loading(in ns)? ex) 200 400 \n")
  tb12=tblock.split()
  tb1,tb2=float(tb12[0]),float(tb12[1])
  ntb1,ntb2=tb1/tstep,tb2/tstep
  ntb1,ntb2=int(round(ntb1)),int(round(ntb2))
  mi12=input("What is the surfactant index interval of MSD calculation? ex) 0 9 \n")
  mi12=mi12.split()
  mi1,mi2=int(mi12[0]),int(mi12[1])

  start_time=timeit.default_timer()

  #input 1 : load surf traj. (big file)
  traj=md.load(trjfile,top=grofile)
  traj=traj[ntb1:ntb2+1]
  topology=traj.topology
  if topology.n_residues <= mi2: #surfactant index exceeded
    mi2=topology.n_residues-1 #autofix surfactant index range
  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))
  print(traj)
  nstep=traj.n_frames

  #prepare bins for MSD
  sdbin,count=numpy.zeros(nstep),numpy.zeros(nstep) #bin & stat weight of dr^2(t) ave

  #input 2 : load migevents. Then, for target surfactant molecules, should make 'cutting time list'.
  # the cutting time list will be used for next loop, to slice trajectories and calc MSD internally for each traj segment.
  #tagmi=mi1 #index of surfactant being searched to write cuttlist
  #row,cuttlist=[],[]
  #while True:
  #  migline=migfile.readline()
  #  if migline=='':#EOF
  #    break
  #  info=mig_stepteller(migline,tb1,tb2,mi1,mi2)
  #  if len(info)!=0 and info[0]==tagmi:
  #    row.append(info[1])
  #  elif len(info)!=0 and info[0]>tagmi:
  #    row.append(nstep)
  #    tagmi=info[0]
  #    cuttlist.append(row)
  #    row=[info[1]]
  #cuttlist.append(row) #for the last molecule\
  total_miginfo=numpy.loadtxt(migfile)
  row,cuttlist=[],[]
  for mi in range(mi1,mi2+1):
    mi_mig=total_miginfo[total_miginfo[:,0]==mi]
    mi_mig=mi_mig[tb1<=mi_mig[:,4]]
    mi_mig=mi_mig[mi_mig[:,4]<=tb2]
    if mi_mig.ndim==1: #1-line
      ntarr=(mi_mig[4]-tb1)/tstep
      ntarr=int(round(ntarr))
      if len(row)==0: row.append(ntarr)
      elif row[-1]!=ntarr: row.append(ntarr) #avoid duplicate
    else:
      for entry in mi_mig:
        ntarr=(entry[4]-tb1)/tstep
        ntarr=int(round(ntarr)) 
        #duplicate check
        if len(row)==0: row.append(ntarr)
        elif row[-1]!=ntarr: row.append(ntarr) #avoid duplicate
    if len(row)==0: row.append(nstep)
    elif row[-1]!=nstep: row.append(nstep)
    if row[0]==0: del row[0]
    cuttlist.append(row)
    row=[]
  
  elapsed=timeit.default_timer() - start_time
  print('migevent information loading complete {}'.format(elapsed))
  print(cuttlist)

  #loop of trajectory slicing and MSD calculation
  mindex=mi1
  for row in cuttlist: #each row represent 1 surfactant
    lastframe=0
    #remove duplicates in 1 row of cuttlist
    for x in row:
      trimtraj=traj[lastframe:x]
      trimtraj=trimtraj.atom_slice(topology.select('resid '+str(mindex)))
      sdbin1,count1=SD_calculator(trimtraj,nstep) #square displacement and statist.count for 1 traj-segment
      sdbin+=sdbin1
      count+=count1
      lastframe=x

      elapsed=timeit.default_timer() - start_time
      print('surf# {} partial trajcut {} MSD calculated. time {:11.4f}'.format(mindex,x,elapsed))
    mindex+=1

  #printing section
  sdbin=numpy.divide(sdbin,count,out=numpy.zeros_like(sdbin), where=count!=0)
  #outfile.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(0,0,0))
  for i in range(nstep):
    outfile.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(tstep*i,sdbin[i],count[i]))

  outfile.close()

if __name__ == "__main__": main()

