#!/home/kjeong23/softwares/bin/python3.4
# umbrella sampling window snapshot detector
# algorithm: read window distance plan file, .xvg or similar format file
# to denote time-distance relation
# load the list of window plan, and use sorting algorithm to find closest distance to each window target

import math
import sys
import numpy
import timeit

def main():
  #Load files
  wipfile = open(sys.argv[1],'r') #window plan
  dstfile = open(sys.argv[2],'r') #dist info file
  outfile = open(sys.argv[3],'w') #output file for corresponding time snapshot
  #env variables
  timecol=0 #column number of time in xvg file
  metric=[4] #list of column number for umbrella sampling coord dis
  start_time=timeit.default_timer()
  #load wipfile -> get info of window design plan
  lindex=0
  for line in wipfile:
    if lindex==0:
      nwin=int(line)
      winlist=numpy.zeros(nwin,float)
    else:
      winlist[lindex-1]=float(line)
    lindex+=1
  #automatic window assignment
  timeinfo=numpy.array([])
  crdinfo=numpy.array([])
  snaptime,actd=numpy.zeros(nwin,float),numpy.zeros(nwin,float)
  for line in dstfile:
    split=line.split()
    if '@' not in split[0] and '#' not in split[0]:
      t=float(split[timecol])
      spsplit=[]
      for x in metric:
        spsplit.append(split[x])
      xvgar=[float(x) for x in spsplit] #raw line segment of xvg file
      xi=numpy.linalg.norm(xvgar)
      timeinfo=numpy.append(timeinfo,t)
      crdinfo=numpy.append(crdinfo,xi)
  for i in range(nwin):
    dev=0.1
    for j in range(len(crdinfo)):
      tempdev=abs(winlist[i]-crdinfo[j])
      if tempdev < dev:
        dev=tempdev
        snaptime[i]=timeinfo[j]
        actd[i]=crdinfo[j]
  elapsed=timeit.default_timer() - start_time
  print('Sorting complete time {}'.format(elapsed))
  #printing section
  outfile.write('  intend    ps    actual    win#\n')
  for i in range(nwin):
    outfile.write('{:8.4f} {:8.4f} {:8.4f} {:3}\n'.format(winlist[i],snaptime[i],actd[i],i))
  wipfile.close()
  dstfile.close()
  outfile.close()
if __name__ == "__main__": main()
