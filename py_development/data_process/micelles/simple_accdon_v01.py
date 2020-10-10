#!/home/kjeong23/softwares/bin/python3.4
# simple lattice site acc-don tracker
# input file : 1 micelle morphology list file 

import sys
import numpy
import timeit

def diff(first, second):
  second = set(second)
  return [item for item in first if item not in second]

def accdon_check(line1,line2):
  memline1,memline2=line1[96:].replace(",",""),line2[96:].replace(",","")
  memsplit1,memsplit2=memline1.split(),memline2.split()
  time1,time2=line1[10:18],line2[10:18]
  accmem=diff(memsplit2,memsplit1)
  donmem=diff(memsplit1,memsplit2)

  info=[time1,time2,len(accmem),len(donmem)]
  return info,accmem,donmem

def main():
  #Load input files
  morfile = open(sys.argv[1],'r')
  outfile = open(sys.argv[2],'w') #output file name.

  #nnmic=int(input("Most frequent number of micelles desired? ex) 30 or 8\n"))
  start_time=timeit.default_timer()

  sindex,acc,don=0,0,0
  oldline=morfile.readline()
  writingcand=[]
  while True:
    newline=morfile.readline()
    if newline=='':#EOF
      break
    info,accmem,donmem=accdon_check(oldline,newline)
    if (info[2]+info[3])!=0:
      row=[info,accmem,donmem]
      writingcand.append(row)
    #at end of step, copy line
    oldline=newline[:]
    sindex+=1 #at the end of 1step loop, add step number

  nwrite=len(writingcand)
  for i in range(nwrite-1): #refine the list of writing. if it's donation/reacception, discard it.
    row1,row2=writingcand[i],writingcand[i+1]
#    time2,time3=float(row1[0][1]),float(row2[0][0])
#    time2,time3=int(time2*1000),int(time3*1000)
#    if time2==time3:
    donmem1,accmem2=row1[2],row2[1]
    common=list(set(donmem1).intersection(accmem2))
    if len(common)!=0:
      writingcand[i][2]=diff(writingcand[i][2],common)
      writingcand[i+1][1]=diff(writingcand[i+1][1],common)
      writingcand[i][0][3]-=len(common)
      writingcand[i+1][0][2]-=len(common)
    accmem1,donmem2=row1[1],row2[2]
    common2=list(set(accmem1).intersection(donmem2))
    if len(common2)!=0:
      writingcand[i][1]=diff(writingcand[i][1],common2)
      writingcand[i+1][2]=diff(writingcand[i+1][2],common2)
      writingcand[i][0][2]-=len(common2)
      writingcand[i+1][0][3]-=len(common2)

  for row in writingcand:
    if row[0][2]+row[0][3]!=0:
      outfile.write('time {} ns to {} ns. accepted {} donated {} '.format(row[0][0],row[0][1],row[0][2],row[0][3]))
      outfile.write(' accmem {} donmem {}\n'.format(row[1],row[2]))
      acc+=int(row[0][2])
      don+=int(row[0][3])

  outfile.write('total acceptance {} donation {} in whole traj for this micelle\n'.format(acc,don))

  morfile.close()
  outfile.close()
   
if __name__ == "__main__": main()
