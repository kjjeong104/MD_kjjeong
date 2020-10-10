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

  aname1=input("atom names for donors? ex) Nm1 Nm2 \n")
  aname2=input("atom names for hydrogens? ex) Hm1 Hm2 Hm3 Hm4 \n")
  aname3=input("atom names for acceptors? ex) CL \n")
  rrange=input("minimum, maximum H-A distance in nanometers? ex) 0.10 0.50 \n")
  rsplit=rrange.split()
  rmin,rmax=float(rsplit[0]),float(rsplit[1])
  rbin=float(input("r bin size? ex) 0.002 \n"))
  abin=float(input("cosine(theta) bin size? ex) 0.02 \n"))
  rnbin,anbin=int((rmax-rmin)/rbin),int(2/abin)
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
  if nstep>=100: #if there're many frames..
    nfrag=200
  else:
    nfrag=1

  elapsed=timeit.default_timer() - start_time
  print('finished trajectory loading {}'.format(elapsed))
  print(nstep,nmon)

  #prepare 2dbins for hydrogen-acceptor distance and H-bond angle
  #sdbin,count=numpy.zeros(nstep),numpy.zeros(nstep) #bin & stat weight of dr^2(t) ave

  #make atom indices list (before filtering too far pairs)
  #should avoid intramolecular atomic pair
  asplit1,asplit2,asplit3=aname1.split(),aname2.split(),aname3.split()
  text1,text2,text3='','',''
  for word in asplit1:
    text1+='name '+word+' or '
  for word in asplit2:
    text2+='name '+word+' or '
  for word in asplit3:
    text3+='name '+word+' or '
  text1,text2,text3=text1[:-4],text2[:-4],text3[:-4]
  seld=topology.select(text1)
  selh=topology.select(text2)
  sela=topology.select(text3)
  n_atomd,n_atomh,n_atoma=len(seld),len(selh),len(sela)
  print(n_atomd,n_atomh,n_atoma)
  dhpairs=[]
  for j in topology.bonds:
    if j[0].index in selh and j[1].index in seld:
      dhpairs.append([j[1].index,j[0].index])
    elif j[0].index in seld and j[1].index in selh:
      dhpairs.append([j[0].index,j[1].index])
  fulllist_angles=[]
  for row in dhpairs:
    for i in sela:
      if topology.atom(row[1]).residue!=topology.atom(i).residue:
        extrow=row.copy()
        extrow.append(i)
        fulllist_angles.append(extrow)
  #list_dist=numpy.array(list_dist)
  fulllist_angles=numpy.array(fulllist_angles)
  n_angles_full = len(fulllist_angles)
  print(" full list # angles = {} ".format(n_angles_full))

  for ifrag in range(nfrag): #loop of fragmental calculation, to save memory.
    blength=int(nstep/nfrag)
    bstart,bend=ifrag*blength,(ifrag+1)*blength
    fragtraj=traj[bstart:bend]
    #calculate distances between hydrogen and acceptors, angle 
    dist = (md.compute_distances(fragtraj,fulllist_angles[:,[1,2]])).flatten()
    #index_thres=numpy.where(full_dist<=rmax)[0]
    #dist=numpy.array(full_dist[index_thres])
    #print(dist)
    #list_angles=fulllist_angles[index_thres,:]
    ## calculate angles and distances
    angl = (md.compute_angles(fragtraj,fulllist_angles)).flatten()
    #n_angles=len(list_angles)
    #print(" list # distances(within threshold) = {}".format(len(dist)))
    #print(" list # angles = {}".format(n_angles))

    # recalculate NaN values of angles
    import copy
    for nan_item in numpy.argwhere(numpy.isnan(angl)).reshape(-1):
        i_frame = int(nan_item/n_angles_full)
        i_angle = nan_item%n_angles_full
        #print(" Nan at {} th frame and {} th angle".format(i_frame,i_angle))
        #print(" <-- {} th atoms".format(list_angles[i_angle]))
        i_abc = fulllist_angles[i_angle]
        a = traj.xyz[i_frame][i_abc[0]]
        b = traj.xyz[i_frame][i_abc[1]]
        c = traj.xyz[i_frame][i_abc[2]]
        print(" <-- position: \n {} \n {} \n {}".format(a,b,c))
        ba = a - b
        bc = c - b
        cosine_angle = numpy.dot(ba, bc) / (numpy.linalg.norm(ba) * numpy.linalg.norm(bc))
        print(" distance= {}".format(numpy.linalg.norm(bc)))
        angle = numpy.arccos(cosine_angle)
        print(" get correct value from NaN, {} (rad) {} (deg)".format(angle, angle*180.0/numpy.pi))
        angl[nan_item] = copy.copy(angle)

    cosangl=numpy.cos(angl)
    #print(cosangl)
    #printing section - should regard gnuplot pm3d-compatible format.
    #hold x. increment y. when a full cycle of y range ends, make an empty line.
    #histogram
    counts_2d, edge_r, edge_cosa = numpy.histogram2d(dist,cosangl,bins=[rnbin,anbin],range=[[rmin,rmax],[-1.0,1.0]])
    #volume in each radial shell
    vol = numpy.power(edge_r[1:],3) - numpy.power(edge_r[:-1],3)
    vol *= 4/3.0 * numpy.pi
    # Average number density
    box_vol = numpy.average(fragtraj.unitcell_volumes)
    density = n_angles_full / box_vol
    rdf_2d = (counts_2d * anbin/ nstep ) / (density* vol[:,None] )
    if ifrag==0:
      totrdf_2d=numpy.copy(rdf_2d)
    else:
      totrdf_2d+=rdf_2d
    elapsed=timeit.default_timer() - start_time
    print('finished fragment {} time {}'.format(ifrag,elapsed))

  for i in range(rnbin):
    for j in range(anbin):
      xval,yval=rmin+rbin*i,-1.0+abin*j
      outfile.write('{:11.4f} {:11.4f} {:11.4f}\n'.format(xval,yval,totrdf_2d[i][j]))
    outfile.write('\n')

  elapsed=timeit.default_timer() - start_time
  print('finished job {}'.format(elapsed))
  #numpy.savetxt(outfile,numpy.transpose(rdf_2d), \
  #header='x = distance [{},{}], y= cos angle [{},{}]' \
  #.format(rmin,rmax,-1.0,1.0),fmt='%f',comments='# ')
  #numpy.save(outfile,numpy.transpose(rdf_2d))
  #print(rdf_2d)
  outfile.close()

if __name__ == "__main__": main()

