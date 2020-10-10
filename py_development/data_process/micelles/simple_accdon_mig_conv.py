#!/home/kjeong23/softwares/bin/python3.4
# acc-don conversion file from migevents
# input file : inter-micellat exchange summary file (migevents.txt)

import sys
import numpy
import timeit

#def diff(first, second):
#  second = set(second)
#  return [item for item in first if item not in second]
#def accdon_check(line1,line2):
#  memline1,memline2=line1[96:].replace(",",""),line2[96:].replace(",","")
#  memsplit1,memsplit2=memline1.split(),memline2.split()
#  time1,time2=line1[10:18],line2[10:18]
#  accmem=diff(memsplit2,memsplit1)
#  donmem=diff(memsplit1,memsplit2)
#  info=[time1,time2,len(accmem),len(donmem)]
#  return info,accmem,donmem

def main():
  #Load input files
  migfile = sys.argv[1]
  outstr = sys.argv[2]

  tstep=0.02
  nmic=30 # number of micelles (preset.)

  start_time=timeit.default_timer()
  #input 1: read migevents. can use numpy array to read everything.
  #while True:
  #  newline=migfile.readline()
  #  if newline=='':#EOF
  #    break
  migevents=numpy.loadtxt(migfile)

  for i in range(1,nmic+1): #loop that write accdon file (nmic) times
    outfile = open(outstr+'_'+str(i)+'.dat','w')
    
    #from total migevents: extract acc/don for corresponding micelle only.
    #then, sort by time.
    mic1_event1=migevents[migevents[:,1]==i]
    mic1_event2=migevents[migevents[:,2]==i]
    mic1_event=numpy.vstack((mic1_event1,mic1_event2))
    #for donation, only take initial time. for acception, only take final time.
    mic1_event_trim=numpy.empty((0,4),float)
    for row in mic1_event: #time donmic accmic surfindex
      if row[1]==i: trimrow=numpy.array([row[3],row[1],row[2],row[0]]) #donation
      elif row[2]==i: trimrow=numpy.array([row[4],row[1],row[2],row[0]])   #acception
      mic1_event_trim=numpy.vstack((mic1_event_trim,trimrow))
    mic1_event_trim=mic1_event_trim[mic1_event_trim[:,0].argsort()]

    #loop writing section
    tacc,tdon=0,0
    for row in mic1_event_trim:
      if row[1]==i: #donation
        t1,t2,nacc,ndon,accmem,donmem=row[0],row[0]+tstep,0,1,"[]",int(row[3])
        tdon+=1
      elif row[2]==i: #acception
        t1,t2,nacc,ndon,accmem,donmem=row[0]-tstep,row[0],1,0,int(row[3]),"[]"
        tacc+=1
      outfile.write('time {:11.4f} ns to {:11.4f} ns. accepted {} donated {}  accmem {} donmem {}\n'.format(t1,t2,nacc,ndon,accmem,donmem))
    outfile.write('total acceptance {} donation {} in whole traj for this micelle\n'.format(tacc,tdon))
    elapsed=timeit.default_timer()-start_time
    print('completed writing file# {}, time {}'.format(i,elapsed))
    outfile.close()

if __name__ == "__main__": main()
