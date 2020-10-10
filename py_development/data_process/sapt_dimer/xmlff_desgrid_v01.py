#!/usr/bin/env python3
# automatic writer for OpenMM .xml force field file, from xml template
# read template, settings
#caution : the 'skeleton' template ff file is the invariable part of total ff file (before customNB definition)
#the last '</ForceField>' line should be missing, because it will be manually printed afterwards.
# the 'rule' file format : atomictype Aexch Aelec Aind Adhf (in atomic units)
# v01: test if writing is correct for a single input.
# v02 plan: grid style workspace generation.

import math
import numpy as np
import argparse
import re

conv_atom_kjmol=2625.5

Desti = {'class' : -1, 'Aexch' : 0, 'Aelec' : 1, 'Aind' : 2, 'Adhf' : 3, 'Bexp' : 4,
	 'C6' : 5, 'C8' : 6, 'C10' : 7, 'C12' : 8, 'Switch' : 9,
	 'Aexadj' : 10, 'Aeladj' : 11, 'Ainadj' : 12, 'Adhadj' : 13}

parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='generate openMM ff xml file (SAPT-FF reline) from template xml FF')
## args
parser.add_argument('-t', '--temp', default='template_ff.xml', nargs='?', 
  help='input rest of xml file template')
parser.add_argument('-i', '--input', default='template_cusnb_ff.xml', nargs='?', 
  help='input customNB partial xml file template')
parser.add_argument('-r', '--rule', default='rule.txt', nargs='?', 
  help='rule file for ff replacement')
parser.add_argument('-o', '--output', default='output_ff.xml', nargs='?', 
  help='output force field file')

def xml_cnb_replace(lines,rule):
  newcnblines = lines.copy()
  data = [line.split() for line in lines] #2d list of the textfile
  for row in rule: #manual replacement check
    rowdata = row.split()
    atype = rowdata[0]
    rnum=0
    #ff file check
    for drow in data:
      if drow[1]=='class="'+atype+'"':
        coefs = np.array([float(x) for x in rowdata[1:5]])*conv_atom_kjmol
        print(coefs)
        coefstr=['{:.1f}'.format(x) for x in coefs]
        #regular expression detection, replacement
        newcnblines[rnum]=re.sub(r"Aexch=\".*\" Aelec","Aexch=\""+coefstr[0]+"\" Aelec",newcnblines[rnum])
        newcnblines[rnum]=re.sub(r"Aelec=\".*\" Aind","Aelec=\""+coefstr[1]+"\" Aind",newcnblines[rnum])
        newcnblines[rnum]=re.sub(r"Aind=\".*\" Adhf","Aind=\""+coefstr[2]+"\" Adhf",newcnblines[rnum])
        newcnblines[rnum]=re.sub(r"Adhf=\".*\" B","Adhf=\""+coefstr[3]+"\" B",newcnblines[rnum])
      rnum+=1

  return newcnblines

# parse xml input file
args = parser.parse_args()

f = open(args.input,'r')
tf = open(args.temp,'r')
rf = open(args.rule,'r')
skellines = tf.readlines()

lines = f.readlines()
rlines = rf.readlines()
newcnblines = xml_cnb_replace(lines,rlines)

#rewrite ff file
outfile=open(args.output,'w')
for x in skellines:
  outfile.write(x)
for x in newcnblines:
  outfile.write(x)
outfile.write('  </CustomNonbondedForce>\n')
outfile.write('</ForceField>\n')

outfile.close()
