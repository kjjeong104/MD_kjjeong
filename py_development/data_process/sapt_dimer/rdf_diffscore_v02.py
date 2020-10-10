#!/usr/bin/env python3
# script for RDF difference score calculation.
# v01: manually put work directories for comparisons.
# algorithm : read xvg file name list with range and inverse weight factor written (gmax or 2*gmax)
# input the reference and comparison work directories.
# rule file format : xvgfilename rcut weight
# v02: tailored for automated grid test. get scan rule
# output rdf score file format: (columns of scan variable) (rdfscore)

import math
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  description='generate gromacs tabulated potential from openmm xml FF')
# args
#parser.add_argument('-i', '--input', default='', nargs='?',
#  help='input atomtypes')
parser.add_argument('-i', '--input', default='', nargs='?',
  help='input directory name pattern.')
parser.add_argument('-n', '--numscan', default='1000', nargs='?',
  help='number of gridscan trials')
parser.add_argument('-f', '--ref', default='', nargs='?',
  help='reference rdf directory name')
parser.add_argument('-r', '--rule', default='rdflist.txt',nargs='?',
  help='rdf names and hyperparameter list')
parser.add_argument('-s', '--scanrule', default='scanrule.txt',nargs='?',
  help='original scanrule for the grid test')
parser.add_argument('-o', '--output', default='output_rdfscore.dat',nargs='?',
  help='output file name for rdf score of the grid test species')
args = parser.parse_args()

#load rdf file
def xvgfile_nparray(fname):
  odata=np.loadtxt(fname,dtype=float,comments=['&','#','@'])
  return odata  

#RDF Mean Square Difference Calculation
def rdf_differ(rdfref,rdf1,rcut): 
  rrange = np.argwhere(rdfref[:,0]<=rcut).flatten()
  #print(rrange)
  prdfref,prdf1 = rdfref[rrange], rdf1[rrange]
  #print(prdfref,prdf1)
  j_rdf = np.average((prdf1[:,1]-prdfref[:,1])**2)
  return j_rdf

# read rdf filename and hyperparameter info
idir_prefix,rdir=args.input,args.ref
ngrid=int(args.numscan)
rf = open(args.rule,'r')
rlines = rf.readlines()
sf = open(args.scanrule,'r')
slines = sf.readlines()
basedir=os.getcwd()
outfile=open(args.output,'w')

ndim=len(slines)
unziprule=np.empty((0,3))
for row in slines:
  rsplit=row.split()
  v1,v2,step=float(rsplit[3]),float(rsplit[4]),float(rsplit[5])
  unziprule=np.vstack((unziprule,np.array([v1,v2,step])))
i=0
if len(slines)==3:
  scanmap=np.empty((ngrid,3))
  for c0 in np.arange(unziprule[0,0],unziprule[0,1],unziprule[0,2]):
    for c1 in np.arange(unziprule[1,0],unziprule[1,1],unziprule[1,2]):
      for c2 in np.arange(unziprule[2,0],unziprule[2,1],unziprule[2,2]):
        scanmap[i]=np.array([c0,c1,c2])
        i+=1

#loop of sequential rdf file comparison
for igrid in range(ngrid):
  idir=idir_prefix+str(igrid)
  total_j_rdf=0.0
  for rulerow in rlines:
    rsplit=rulerow.split()
    fname,rcut,invweight=rsplit[0],float(rsplit[1]),float(rsplit[2])
    desti_f,desti_i = os.path.join(basedir,rdir,fname),os.path.join(basedir,idir,fname)
    rdfref_data,rdf1_data=xvgfile_nparray(desti_f),xvgfile_nparray(desti_i)
    j1_rdf = rdf_differ(rdfref_data,rdf1_data,rcut)
    total_j_rdf += j1_rdf/invweight

  #print (rdir,idir,total_j_rdf)
  outfile.write('{:.3f} {:.3f} {:.3f} {:9.4f}\n'.format(scanmap[igrid,0],scanmap[igrid,1],scanmap[igrid,2],total_j_rdf))
  if igrid%50==0:
    print('evaluation complete: '+str(igrid))

outfile.close()

