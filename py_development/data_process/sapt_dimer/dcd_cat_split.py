#!/home/kjeong23/softwares/bin/python3.4
# program to concatenate dcd trajectories, split it if needed
#inputs: pdbfile, traj dcd file(s).
#output: !!Do not put output file name as argument !! fragment dcd file

import math
import sys
import numpy
import mdtraj as md

#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  alength = len(sys.argv) #number of dcd files can vary
  #trjfile = sys.argv[2]
  #outfile = sys.argv[3]

  #settings for env variables.
  outstr=input("output dcd file prefix? ex) DES_npt_b (becomes DES_npt_b${i}.dcd) \n")
  nb=int(input("how many blocks? ex) 5 \n"))
  bstep=int(input("frames per block? ex) 2000 \n"))
  neq=int(input("How many initial frames do you want to cut as equilibration? ex) 500 \n"))

  #input 1 : load traj. (big file)
  #check if multiple trajectories are put and needed to be merged.
  ntrajfile=alength-2
  if ntrajfile==1: 
    trajtotal=md.load(sys.argv[2],top=topfile)
    nframetotal=trajtotal.n_frames
    print('nframetotal : {}'.format(nframetotal))
  elif ntrajfile==2:
    traj1,traj2=md.load(sys.argv[2],top=topfile),md.load(sys.argv[3],top=topfile)
    nframe1,nframe2=traj1.n_frames,traj2.n_frames
    print('nframe1, nframe2 : {} {}'.format(nframe1,nframe2))
    trajtotal=md.join([traj1,traj2],check_topology=True,discard_overlapping_frames=True)
    nframetotal=trajtotal.n_frames
    print('nframetotal : {}'.format(nframetotal))

  #traj splitting:loop
  for i in range(nb):
    outfile=outstr+str(i)+'.dcd'
    iframe,fframe=neq+i*bstep,neq+(i+1)*bstep
    trajfrag=trajtotal[iframe:fframe]
    trajfrag.save_dcd(outfile)

if __name__ == "__main__": main()

