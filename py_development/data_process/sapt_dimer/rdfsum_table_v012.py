#!/usr/bin/env python3
#python script file for comprehensive rdf summary table generation
#columns : rdftype rmax1 rmin1 peakheight  CN
# algorithm : set loop of mutual atom types. in the atomtype loop:
# read xvg file of gr, detect rmax1, rmin1, peakheight. 
# read cn file. according to rmin, read CN of 1st solvation shell.
# write all values in an output table.

import math
import numpy as np
from scipy.signal import find_peaks
import sys
import os.path
import argparse
parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='generate gromacs tabulated potential from openmm xml FF')
# args
#parser.add_argument('-i', '--input', default='', nargs='?', 
#  help='input atomtypes')
parser.add_argument('-o', '--output', default='rdfsum_table.txt', nargs='?', 
  help='input customNB partial xml file')
parser.add_argument('-w', '--width', default=5,nargs='?',
  help='peak width criterion for local extrema detection')
args = parser.parse_args()

def xvgfile_nparray(fname):
  skipnum=0
  with open(fname) as f:
    for line in f:
      if '#' in line or '@' in line:
        skipnum+=1
      else:
        break
  #print('debug : skipnum = ',skipnum)
  odata=np.loadtxt(fname,dtype=float,skiprows=skipnum)
  return odata

def hydshell_report(grdata,cndata,w):
  w=int(w)
  peaks,_ = find_peaks(grdata[:,1],distance=w,width=w)
  revpeaks,_ = find_peaks(-grdata[:,1],distance=w,width=w)
  
  imax1,imin1=peaks[0],revpeaks[0]
  rmax1,rmin1,peakheight=grdata[imax1][0],grdata[imin1][0],grdata[imax1][1]
  cn=cndata[imin1][1]   
 
  hydshell_line=[rmax1,rmin1,peakheight,cn]
  return hydshell_line
######################################################

outfile = open(args.output,'w')
w=args.width
outfile.write('{:10} {:10} {:9} {:9} {:9} {:9}\n'.format('atom1','atom2','rmax','rmin','peak','CN'))
atypestr = input("List all atom types to inspect. ex) Cl uH cholOH \n")
atypes = atypestr.split()
fstrin = input("Filename prefix without underbar for gr,cn xvg files. ex) gr cn \n")
prefstr = fstrin.split()
for at1 in atypes:
  for at2 in atypes:
    fstr1=prefstr[0]+'_'+at1+'_'+at2+'.xvg'
    if os.path.isfile(fstr1): #grfile exists
      fstr2=prefstr[1]+'_'+at1+'_'+at2+'.xvg'
      #grfile reading, numpy processing. need to sort out string lines.
      grdata=xvgfile_nparray(fstr1)
      cndata=xvgfile_nparray(fstr2)
      hydshell_line=hydshell_report(grdata,cndata,w)
      #print(grdata[grmax,:],grdata[grmin,:])
      outfile.write('{:10} {:10} {:9.4f} {:9.4f} {:9.4f} {:9.4f}\n'.format(at1,at2,hydshell_line[0],hydshell_line[1],\
hydshell_line[2],hydshell_line[3]))
    else:
      continue

outfile.close()

