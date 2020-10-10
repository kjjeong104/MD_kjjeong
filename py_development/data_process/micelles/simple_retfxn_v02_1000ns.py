#!/home/kjeong23/softwares/bin/python3.4
# simple time-dependent retention fraction calculation
# input file : 1 micelle acc/don event list file. (made by simple_accdon_v01.py), 1micelle member history file (track)
# Be careful: even if surfactant molecule comes back by acc/don,
# it doesn't count as retention anymore.
# Algorithm: read all acception,donation events and its time, then
# construct a 'history row vector' of aggN,#don with respect to time.
# To calculate retention function, can just 'sum up' tau interval in donation history vector.
# v02 : solve the problem that repetitive donation of same surfactant is overcounted
# modified algorithm : make "donlist" numpy array of don event. (time,surf#). So calc repetition

import sys
import numpy
import timeit

tstep=0.02
nstep=50000

#info : sindex, nacc, ndon.  donevent: row or array that contains 'sindex','donation surfindex'
def split_accdon(line):
  split=line.split()
  time,nacc,ndon=float(split[4]),int(split[7]),int(split[9])
  sindex=time/tstep
  sindex=int(round(sindex))
  info=numpy.array([sindex,nacc,ndon]) 

  subline=line[line.index('donmem')+len('donmem'):]
  subline=subline.replace("[","")
  subline=subline.replace("]","")
  subline=subline.replace("'","")
  subline=subline.replace(",","")
  donsplit=subline.split()

  donevent=numpy.empty((0,2),int)
  for x in donsplit:
    numrow=numpy.array([sindex,int(x)])
    donevent=numpy.vstack((donevent,numrow))
  return info,donevent

def true_donnum(t1,t2,donlist): #calculate 'true' total donation number during a time interval, decline overcount
  #t1 : beginning of interval. Not including itself. t2 : end of interval, include itself.
  time_vector=donlist[:,0]
  member_vector=donlist[:,1]
  i1,i2=numpy.searchsorted(time_vector,t1,side='right'),numpy.searchsorted(time_vector,t2,side='right')
  member_vector=member_vector[i1:i2] #donated molecules indices list with overcounting.
  member_vector=numpy.unique(member_vector)#sort out overcounted
  donnum=len(member_vector)
  return donnum

def main():
  #Load input files
  evefile = open(sys.argv[1],'r') #acc/don event file
  trafile = open(sys.argv[2],'r') #1 micelle tracking file, to sort out unassigned snapshot
  outfile = open(sys.argv[3],'w') #output file name.
  #nsurf=int(input("Initial number of surfactants in this micelle? ex) 32\n"))
  #tstep=0.02
  #nstep=50000
  taumin=float(input("minimum correlation time you want (in ns)? ex) 0\n"))
  taumax=float(input("maximum correlation time you want (in ns)? ex) 100\n"))
  nmintau=int(taumin/tstep)
  nmaxtau=int(taumax/tstep)
  nsteptau=nmaxtau-nmintau

  start_time=timeit.default_timer()

  #arrays
  aggn,acc,don,rf=numpy.zeros(nstep+1),numpy.zeros(nstep+1),numpy.zeros(nstep+1),numpy.zeros(nsteptau)
  donlist=numpy.empty((0,2),int) #donlist: N*2 array of 'sindex,donated surf index'

  #read tracking file, to determine initial aggregation number and mark on unassigned snapshots
  traline=trafile.readline()
  split=traline.split()
  nsurf=int(split[5]) #initial aggN
  aggn[0]=nsurf

  #read eventfile, to write acc and don first. Then, can construct aggn history.
  #(don't use aggn history of trafile directly, because it contains occasional detection interruption) 
  while True:
    newline=evefile.readline()
    if newline[0:5]=='total':#EOF
      break
    info,donevent=split_accdon(newline)
    sindex=info[0]
    acc[sindex],don[sindex]=info[1],info[2]
    if donevent.shape[0]!=0:
      for row in donevent:
        donlist=numpy.vstack((donlist,row))
  #aggn history
  for i in range(1,nstep+1):
    aggn[i]=aggn[i-1]+acc[i]-don[i]

  #read rest of the tracking history file, to delete unassigned snapshots.
  preindex=0
  while True:
    traline=trafile.readline()
    if traline=='':#EOF
      break
    split=traline.split()
    sindex=float(split[2])/tstep
    sindex=int(round(sindex))
    if sindex-preindex>1: #there's unassigned snapshot in between recorded snapshots
      for i in range(preindex+1,sindex):
        aggn[i]=-100 #dummy number
    preindex=sindex

  elapsed=timeit.default_timer() - start_time
  print('Reading input information complete, time {:11.4f}'.format(elapsed))

  #calc retention function: tau loop, initial time loop
  taindex=0
  for j in range(nmintau+1,nmaxtau+1):
    count=0
    for i in range(nstep+1-j):
      if aggn[i]!=-100 and aggn[i+j]!=-100:
        oriagg=aggn[i]
        donnum=true_donnum(i,i+j,donlist)
        #newagg=oriagg-numpy.sum(don[i+1:i+j+1])
        newagg=oriagg-donnum
        if newagg<=0: newagg=0
        frac=float(newagg/oriagg)
        rf[taindex]+=frac
        count+=1
    rf[taindex]/=count #normalize
    taindex+=1
    if j%100==0:
      elapsed=timeit.default_timer() - start_time
      print('time correlating tau-step {}, time {:11.4f}'.format(j,elapsed))

  #output printing section
  if nmintau==0:
    outfile.write('{:9.4f} {:11.6f}\n'.format(0,1))
  for i in range(rf.size):
    tau=taumin+(i+1)*tstep
    outfile.write('{:9.4f} {:11.6f}\n'.format(tau,rf[i]))

  evefile.close()
  trafile.close()
  outfile.close()
   
if __name__ == "__main__": main()
