#!/usr/bin/env python3
# script for RDF difference score calculation.
# v01: manually put work directories for comparisons.
# algorithm : read xvg file name list with range and inverse weight factor written (gmax or 2*gmax)
# input the reference and comparison work directories.
# rule file format : xvgfilename rcut weight

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
  help='input directory name')
parser.add_argument('-f', '--ref', default='', nargs='?',
  help='reference rdf directory name')
parser.add_argument('-r', '--rule', default='rdflist.txt',nargs='?',
  help='rdf names and hyperparameter list')
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
idir,rdir=args.input,args.ref
rf = open(args.rule,'r')
rlines = rf.readlines()
basedir=os.getcwd()

#loop of sequential rdf file comparison
total_j_rdf=0.0
for rulerow in rlines:
  rsplit=rulerow.split()
  fname,rcut,invweight=rsplit[0],float(rsplit[1]),float(rsplit[2])
  desti_f,desti_i = os.path.join(basedir,rdir,fname),os.path.join(basedir,idir,fname)
  rdfref_data,rdf1_data=xvgfile_nparray(desti_f),xvgfile_nparray(desti_i)
  j1_rdf = rdf_differ(rdfref_data,rdf1_data,rcut)
  total_j_rdf += j1_rdf/invweight

print (rdir,idir,total_j_rdf)
