#!/usr/bin/env python3
# automatic writer for OpenMM .xml force field file, from xml template
# read template, settings
#caution : the 'skeleton' template ff file is the invariable part of total ff file (before customNB definition)
#the last '</ForceField>' line should be missing, because it will be manually printed afterwards.
# the 'rule' file format : atomictype Aexch Aelec Aind Adhf (in atomic units)
# v01: test if writing is correct for a single input.
# v02: grid style workspace generation. Also, change the format of ff generation order -> should accept grid scan.
# also, v02 should recognize switchable Cl as 'swCl' atomtype to change 'A adj' parameters.

import math
import numpy as np
import argparse
import re
import os

conv_atom_kjmol=2625.5

Desti = {'class' : -1, 'Aexch' : 0, 'Aelec' : 1, 'Aind' : 2, 'Adhf' : 3, 'Bexp' : 4,
	 'C6' : 5, 'C8' : 6, 'C10' : 7, 'C12' : 8, 'Switch' : 9,
	 'Aexadj' : 10, 'Aeladj' : 11, 'Ainadj' : 12, 'Adhadj' : 13}
Desti2 ={"Aexch": r"Aexch=\".*\" Aelec", "Aelec": r"Aelec=\".*\" Aind", "Aexadj" : r"Aexadj=\".*\" Aeladj"}
Desti3 ={"Aexch": "Aelec", "Aelec": "Aind", "Aexadj": "Aeladj"}

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

def xml_cnb_replace2(lines,rule):
  newcnblines = lines.copy()
  data = [line.split() for line in lines] #2d list of the textfile
  for row in rule: #manual replacement check
    #rowdata = row
    atype,ctype,v = row[0],row[1],row[2]
    rnum=0
    #ff file check
    for drow in data:
      if drow[1]=='class="'+atype+'"':
        coefs = float(v)*conv_atom_kjmol 
        coefstr='{:.1f}'.format(coefs)
        #regular expression detection, replacement
        newcnblines[rnum]=re.sub(Desti2[ctype],ctype+"=\""+coefstr+"\" "+Desti3[ctype],newcnblines[rnum])
        #newcnblines[rnum]=re.sub(r"Aexch=\".*\" Aelec","Aexch=\""+coefstr[0]+"\" Aelec",newcnblines[rnum])
        #newcnblines[rnum]=re.sub(r"Aelec=\".*\" Aind","Aelec=\""+coefstr[1]+"\" Aind",newcnblines[rnum])
        #newcnblines[rnum]=re.sub(r"Aind=\".*\" Adhf","Aind=\""+coefstr[2]+"\" Adhf",newcnblines[rnum])
        #newcnblines[rnum]=re.sub(r"Adhf=\".*\" B","Adhf=\""+coefstr[3]+"\" B",newcnblines[rnum])
      rnum+=1

  return newcnblines

#function to unzip the scan information to list of rules.
def scan_unzip(scanrule):
  unzipped_rule=[]
  atypes,ctypes=[],[]
  coefsave=[]
  rnum=0
  nclrow,nswclrow=-1,-1
  for row in scanrule:
    rsplit=row.split()
    atype,ctype=rsplit[0],rsplit[1]
    if atype=='Cl':
      nclrow=rnum
    elif atype=='swCl':
      nswclrow=rnum
    atypes.append(atype)
    ctypes.append(ctype)
    if rsplit[2]=='scan':
      v1,v2,step=float(rsplit[3]),float(rsplit[4]),float(rsplit[5])
      values=np.arange(v1,v2,step)
      coefsave.append(values)
    rnum+=1
  #switchable Cl exists
  #Cl(pairHo,swCl) = Cl(urea,Cl)+adjCl
  #adjCl = swCl-Cl(urea,Cl)
  if nclrow!=-1 and nswclrow!=-1:
    atypes[nswclrow]='Cl'
    ctypes[nswclrow]=ctypes[nswclrow][0:3]+'adj'
  #ad-hoc solution, if scan variables are exactly 3 dimensions.
  if len(scanrule)==3:
    for c0 in coefsave[0]:
      for c1 in coefsave[1]:
        for c2 in coefsave[2]:
          rule1=[]
          c=np.array([c0,c1,c2])
          #reprocess swCl
          if nclrow!=-1 and nswclrow!=-1:
            c[nswclrow]=c[nswclrow]-c[nclrow]
          rule1.append([atypes[0],ctypes[0],c[0]])
          rule1.append([atypes[1],ctypes[1],c[1]])
          rule1.append([atypes[2],ctypes[2],c[2]])
          unzipped_rule.append(rule1)
  return unzipped_rule

# parse xml input file
args = parser.parse_args()

f = open(args.input,'r')
tf = open(args.temp,'r')
rf = open(args.rule,'r')
skellines = tf.readlines()
lines = f.readlines()
rlines = rf.readlines()

#unzip the scanrule
unzipped_rule=scan_unzip(rlines)

#grid workspace creation, and ff file writing.
basedir=os.getcwd()
dirstr="grid"
i=0
for rule1 in unzipped_rule: 
  newcnblines = xml_cnb_replace2(lines,rule1)
  
  #workspace, print
  desti=os.path.join(basedir,dirstr+str(i))
  if os.path.exists(desti):
    print(desti + ' : exists')
  else:
    os.mkdir(desti)
  os.system('cp *run*.py '+desti)
  os.chdir(desti)
  outfile=open(args.output,'w')
  for x in skellines:
    outfile.write(x)
  for x in newcnblines:
    outfile.write(x)
  outfile.write('  </CustomNonbondedForce>\n')
  outfile.write('</ForceField>\n')
  outfile.close() 

  i+=1

