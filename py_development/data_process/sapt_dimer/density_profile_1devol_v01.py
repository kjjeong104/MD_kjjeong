#!/usr/bin/env python 
# density profile calculator. 1-dimensional projection, time evolution.
# density is expressed in number of molecules per nm^3.
#output data type : 2d contour type. x axis is the primary axis coordinate, y axis is time.
# z axis value is the corresponding number density.

#algorithm : read trajectory and settings(nbin, primary axis, time range, timestep),
# -> detect species. -> COM reduction -> binning the average number density into 1d axis.

#v01 : not compatible with polymers. 1d only supported.
# for npt simulation : approximation : use the minimum box size for the density profile calculation.(avoid 0 entry)
import mdtraj as md
import numpy as np
import timeit
import argparse
import sys

parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='density profile calculator')
# args
parser.add_argument('-s', '--topfile', default='', nargs='?', 
  help='topology pdb file')
parser.add_argument('-f', '--trjfile', default='', nargs='?', 
  help='trajectory file')
parser.add_argument('-o', '--outfile', default='den1dvol', nargs='?', 
  help='dipole moment histogram for the trajectory')
parser.add_argument('-b', '--begin',default=0,nargs='?',
  help='initial cut frames for equilibration ')
parser.add_argument('-skip', '--tskip',default=1,nargs='?',
  help='actual timestep regarding skipping')
parser.add_argument('-n', '--nbin',default=100,nargs='?',
  help='number of spatial bins in density profile')
parser.add_argument('-a', '--axis',default='z',nargs='?',
  help='primary axis')
parser.add_argument('-i', '--interactive',default='no',nargs='?',
  help='interactive additional options')
args = parser.parse_args()

def pm3d_print(outfile,array,xmin,xinc,ymin,yinc): #void function for gnuplot pm3d compatible printing 
  #assume array is 2-dimensional. z values are stored in [nx,ny] size array.
  xnbin,ynbin=array.shape[0],array.shape[1]
  outfilef=open(outfile,'w')
  for i in range(xnbin):
    xval=xmin+i*xinc
    for j in range(ynbin):
      yval=ymin+j*yinc
      outfilef.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(xval,yval,array[i][j]))
    outfilef.write('\n')
  outfilef.close()

def main():
  #Part to load coordinate file
  topfile = args.topfile
  trjfile = args.trjfile
  #urefile = args.urefile
  outstr = args.outfile
  teq,tskip = int(args.begin),int(args.tskip)
  nbin=int(args.nbin)
  if args.axis=='x':
    axisi=0
  elif args.axis=='y':
    axisi=1
  else:
    axisi=2
  #input 1 : load surf traj. (big file)
  #t=md.load(topfile)
  traj=md.load(trjfile,top=topfile)
  topology = traj.topology
  traj=traj[teq::tskip]
  nstep=traj.n_frames
  print(nstep, ' frames ')

  #box information, density bin processing. 
  full_boxinfo=traj.unitcell_lengths
  min_boxinfo=full_boxinfo[ np.argmin(full_boxinfo[:,axisi]) ]
  print(min_boxinfo)
  zmin,zmax=0.0,min_boxinfo[axisi]
  dz = float(zmax-zmin)/nbin
  if axisi==2:
    xsize,ysize=min_boxinfo[0],min_boxinfo[1]
  
  #If interactive mode is on, do additional survey. If a specific atom is required instead of COM.
  iresp = args.interactive
  if iresp == 'yes':
    atlist=input("Any specific atomtype target instead of molecular COM? ex) C2a OW \n").split()

  #biggest loop : species loop
  for atype in atlist:
    traj_species=traj.atom_slice(topology.select('name == '+atype))
    denprof_t1d=np.empty((nstep,nbin)) #data array.
    #warning : axis0 is time, axis1 is spatial here, but it has to be switched.
    #timestep reading loop.
    for i in range(nstep):
      traj_1step=traj_species[i]
      atom_crd_axis=traj_1step.xyz[:,:,axisi].flatten()
      #print(atom_crd_axis)
      hist,bin_edges=np.histogram(atom_crd_axis,bins=nbin,range=[zmin,zmax],density=False)
      #average number density = (number in a slab)/(slab volume) = (num)/(xsize*ysize*dz)
      denprof_t1d[i]=(hist/(xsize*ysize*dz)).copy()

    #output processing
    outfile=outstr+'_'+atype+'.dat'
    pm3d_print(outfile,denprof_t1d.transpose(),zmin,dz,0.0,traj.timestep)

if __name__ == "__main__": main()

