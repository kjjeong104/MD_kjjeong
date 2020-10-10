#!/usr/bin/env python
# an intermediate output processing program from spatd_hbond_counter_v02 code
# reads hydrogen bond 'matrix' file, do post analyses
# matrix file format: rows: snapshots. column: by molecule indices ( record the number of hydrogen bonds at that snapshot )

# algorithm : read dcd file to get coord. also, get donor-hydrogen-acceptor atomtype series.
# ex) HBD molecule: urea (NH), HBA: Cl. then Nm Hm Cl
# determine 'upper bound distance' for hydrogen-acceptor distance (ex 5.0A) to check.
# also determine bin widths for r, cos(theta).
# collect D-A distances within threshold, then make angle calculation atom index list.
# calc cosine of dha angle. collect into 2-dimensional bins
#inputs: pdbfile, traj dcd file. interactive input for settings parameters
#not supporting 'spherical or ecliptical hydrogen bond cutoff definition' yet.
#also, use angstrom for basic length unit (be careful that mdtraj use nanometers as basic distance unit)
# counter_v02 update : from experience, solid identification needs short-time average of molecular n_HB.
# therefore, v02 calculates prob-dis of short-time n_HB average, 
# then calculate number of solid molecules at given time (judging vicinity time snapshots) (trajectory fragmental calculation algorithm revived.)
# ex) Jun Soo Kim's ice research : used 20ps interval average of n_HB, when 0.4ps was raw snapshot timestep. (50 snapshot ave) here try 20ps too.
# in this code: supports calculations with many different short-time interval length for nHB averaging, also nsol cutoff
#v011 : changed according to hbond matrix v011 update. considered both donor and acceptor.
#thbave: do not consider seriously anymore. just use 50 ps. histogram range: from 0~4, now 0~8.(bin 0.2)
#nhbscut_list: now change cutoff inspection: 1.0~4.0 with 0.2 interval

import math
import sys
import numpy
import copy
import argparse
import timeit

parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='hydrogen bond solid identifier')

# args
parser.add_argument('-i', '--inputfile', default='hb_matrix.dat', nargs='?', 
  help='hydrogen bond matrix file')
parser.add_argument('-hs', '--histofile', default='histnhb', nargs='?', 
  help='hydrogen bond matrix file')
parser.add_argument('-s', '--solidfile', default='nsol', nargs='?', 
  help='hydrogen bond matrix file')
parser.add_argument('-t', '--tstep', default='2.0', nargs='?', 
  help='snapshot timestep in picoseconds')
args = parser.parse_args()

#main fxn
def main():
  hbmatfile = args.inputfile
  nhboutfile = args.histofile #short time n_HB ave outfile
  nsoutfile = args.solidfile  #number of solid molecules outfile (time evolution)
  tstep=float(args.tstep)

  thbave_list=numpy.array([x*10.0 for x in range(5,6)]) #short-time n_HB average candidates
  nhbscut_list=numpy.array([x*0.2 for x in range(5,21)]) #Desired short-time n_HB cutoff for solid identification
  #thbave=float(input("How many ps to calculate short-time n_HB average? ex) 20 \n"))
  #nhbscut=float(input("Desired short-time n_HB cutoff for solid identification? ex) 2.50 \n")) 

  #reload hydrogen bond matrix. (note that hbmat file is in transposed form. revert it back.)
  hb_total_matrix=numpy.loadtxt(hbmatfile)
  mon_hbond_count=numpy.transpose(hb_total_matrix)
  nmon,nstep=mon_hbond_count.shape[0],mon_hbond_count.shape[1]
  print(nstep,'steps',nmon,'monomers')

  for thbave in thbave_list: #determined thbave first.
    nthbave=int(thbave/tstep)
    nstint=int(nstep/nthbave)
    nhboutfile_finalname=nhboutfile+'_stps_'+str(int(thbave))+'.dat'
    print(nthbave,' number of snapshots for nHB averaging')
    stave_nhb=numpy.empty(nstint*nmon)
    stint=numpy.zeros(nstint)
    for k in range(nstint):
      stint[k]=thbave*k
      for imon in range(nmon):
        stave=numpy.average(mon_hbond_count[imon][k*nthbave:(k+1)*nthbave])
        stave_nhb[k*nmon+imon]=stave
    #histogram drawing
    nbin=40
    if nthbave==5:
      nbin=20
    hist_stave_nhb,bin_edges=numpy.histogram(stave_nhb,bins=nbin,range=[0.0,8.0],density=True)
    bin_edges=bin_edges[:-1]
    bin_edges,hist_stave_nhb= bin_edges.reshape(-1,1),hist_stave_nhb.reshape(-1,1)
    hist_stave_nhb=numpy.hstack((bin_edges,hist_stave_nhb))
    numpy.savetxt(nhboutfile_finalname,hist_stave_nhb,fmt='%.3e')

    #solid identification cutoff test
    for nhbscut in nhbscut_list:
      nsoutfile_finalname=nsoutfile+'_stps_'+str(int(thbave))+'_scut_'+'{:3.1f}'.format(nhbscut)+'.dat'
      nsol=numpy.zeros(nstint)
      for k in range(nstint):
        target=stave_nhb[k*nmon:(k+1)*nmon]
        nsol[k]+=numpy.where(target >=nhbscut)[0].shape[0]
      nsol_display=numpy.hstack((stint.reshape(-1,1),nsol.reshape(-1,1)))
      numpy.savetxt(nsoutfile_finalname,nsol_display,fmt='%d')

if __name__ == "__main__": main()

