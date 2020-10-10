#!/home/kjeong23/softwares/bin/python3.4
# program to calculate distance-angle distribution function of hbond trio.
# algorithm : read dcd file to get coord. also, get donor-hydrogen-acceptor atomtype series.
# ex) HBD molecule: urea (NH), HBA: Cl. then Nm Hm Cl
# determine 'upper bound distance' for hydrogen-acceptor distance (ex 5.0A) to check.
# also determine bin widths for r, cos(theta).
# collect H-A distances within threshold, then make angle calculation atom index list.
# calc cosine of dha angle.
# collect into 2-dimensional bins
#inputs: pdbfile, traj dcd file
#output: file of {r(H-A), cos(theta), probability density}
# be careful that mdtraj use nanometers as basic distance unit

import math
import sys
import numpy
import timeit
import mdtraj as md

#tstep=1.000 #used ps time unit. 1 snapshot : 1ps.
#main fxn
def main():
  #Part to load coordinate file
  topfile = sys.argv[1]
  trjfile = sys.argv[2]
  outfile = open(sys.argv[3],'w')

  aname1=input("atom names for group 1 ex) Hm1 Hm2 Hm3 Hm4 \n")
  aname2=input("atom names for group 2 ex) CL \n")
  rrange=input("minimum, maximum H-A distance in nanometers? ex) 0.10 0.50 \n")
  rsplit=rrange.split()
  rmin,rmax=float(rsplit[0]),float(rsplit[1])
  rbin=float(input("r bin size? ex) 0.002 \n"))
  rnbin=int((rmax-rmin)/rbin)
  #mi12=input("What is the residue number index interval of MSD calculation? ex) 0 299 \n")
  #mi12=mi12.split()
  #mi1,mi2=int(mi12[0]),int(mi12[1])
  tskip=int(input("Once in how many frames do you want to take? ex) 10 \n"))
  teq=int(input("How many initial frames do you want to cut as equilibration? ex) 5000 \n"))

  start_time=timeit.default_timer()

  #input 1 : load surf traj. (big file)
  traj=md.load(trjfile,top=topfile)
  traj=traj[teq::tskip]
  topology=traj.topology
  nstep=traj.n_frames
  nmon=topology.n_residues

  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))
  print(nstep,nmon)

  #prepare 2dbins for hydrogen-acceptor distance and H-bond angle
  #sdbin,count=numpy.zeros(nstep),numpy.zeros(nstep) #bin & stat weight of dr^2(t) ave

  #make atom indices list (before filtering too far pairs)
  asplit1,asplit2=aname1.split(),aname2.split()
  text1,text2='',''
  for word in asplit1:
    text1+='name '+word+' or '
  for word in asplit2:
    text2+='name '+word+' or '
  text1,text2=text1[:-4],text2[:-4]
  sel1=topology.select(text1)
  sel2=topology.select(text2)
  n_atom1,n_atom2=len(sel1),len(sel2)
  print(n_atom1,n_atom2)
  fulllist_dist=[]
  for i in sel1:
    for j in sel2:
      fulllist_dist.append([i,j])
      
  #list_dist=numpy.array(list_dist)
  fulllist_dist=numpy.array(fulllist_dist)
  n_dist_full = len(fulllist_dist)
  print(" full list # dist = {} ".format(n_dist_full))

  #calculate distances between hydrogen and acceptors, angle 
  full_dist = (md.compute_distances(traj,fulllist_dist)).flatten()

  counts,edge_r = numpy.histogram(full_dist,bins=rnbin,range=[rmin,rmax])
  #volume in each radial shell
  vol = numpy.power(edge_r[1:],3) - numpy.power(edge_r[:-1],3)
  vol *= 4/3.0 * numpy.pi
  # Average number density
  box_vol = numpy.average(traj.unitcell_volumes)
  density = n_dist_full / box_vol
  rdf = (counts / nstep ) / (density* vol)

  for i in range(rnbin):
    xval=rmin+rbin*i
    outfile.write('{:11.4f} {:11.4f}\n'.format(xval,rdf[i]))

  #for i in range(rnbin):
  #  for j in range(anbin):
  #    xval,yval=rmin+rbin*i,-1.0+abin*j
  #    outfile.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(xval,yval,rdf_2d[i][j]))
  #  outfile.write('\n')
  outfile.close()

if __name__ == "__main__": main()

