#!/usr/bin/env python
# OpenMM force field xml file modifier. scale polarizability, by percent of alpha value.
# reduced formatting of xml file : reading entire file is hard.
# parts needed to be changed due to polarizability change ; nonbondedforce charge value, drudeforce charge value, polarizability value.
#syntax : feed (recognize 'D' entries) NB section and Drude section of full-pol.
#algorithm : calculate proportionality const in alpha = q^2/k. calculate parent atom charge. 
#rescale alpha, rescale q, then determine all charges.
#in the input file : skip drudeless atoms.

import numpy as np
import argparse

parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='modify polarizability from partially written openmm xml FF')
## args
parser.add_argument('-i', '--input', default='ff.xml', nargs='?', 
  help='input customNB partial xml file')
parser.add_argument('-o', '--output', default='ffout.xml', nargs='?', 
  help='output customNB partial xml file')
parser.add_argument('-p','--percent',default='0.0',nargs='?',
  help='change of polarizability in percent')

# parse xml input file
args = parser.parse_args()

f = open(args.input,'r')
lines = f.readlines()
data = [line.split() for line in lines]
natoms = len(data)
coefs=np.zeros((natoms,14))
aclass=[None]*natoms


# 

#construct parent atom net charge array

#construct drude polarizability array

#construct drude charge array

#calculate k for alpha= q^2/k.

#construct new alpha, new q.

#redistribute parent atom charge, with writing drude charge

for i in range(natoms):
  line=data[i]
  for word in line:
    if '=' in word:
      wsplit=word.split('=')
      prefix=wsplit[0]
      value=wsplit[1].split('"')[1]
      k=Desti[prefix]
      if k==-1:
        aclass[i]=value
      else:
        coefs[i][k]=np.float64(value)
 
