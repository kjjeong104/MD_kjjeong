#!/usr/bin/env python
#inputs: input ASCII file, row replacement rule file (2 columns: past rownum, new rownum)
#parameters: initial shift linenum. periodicity.
#output: processed ASCII file

import sys
import numpy as np
from itertools import islice

#main fxn
def mapping(lines,rule): #lines: lines of text, rule: mapping rule array
  maplines=lines.copy()
  for onemap in rule:
    i0,i1=int(onemap[0]-1),int(onemap[1]-1)
    maplines[i1]=lines[i0]
    #maplines[i0]=lines[i1]
  return maplines

def main():
  #Part to load coordinate file
  inpfile=sys.argv[1]
  rulefile=sys.argv[2]
  outfile=open(sys.argv[3],'w')
  rule=np.loadtxt(rulefile)

  shift=int(input("Initial shift line number? ex) 2 \n"))
  period=int(input("periodicity? \n"))
 
  #big loop of whole inputfile reading
  with open(inpfile) as f:
    init_nlines = list(islice(f,shift))
    for row in init_nlines:
      outfile.write(row)
    while True:
      lines=list(islice(f,period))
      if not lines:
        break
      maplines=mapping(lines,rule)
      for row in maplines:
        outfile.write(row)

if __name__ == "__main__": main()

