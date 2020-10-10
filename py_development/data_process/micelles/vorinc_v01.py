#!/home/kjeong23/softwares/bin/python3.4
# Voronoi polyhedron counterion/water inclusion counting program
# algorithm: Don't actually calculate Voronoi polyhedron! 
# get atomindex info from grofile -> read COMfile(1step) and ASCII traj(1step) simultaneously
# -> classify each MeN, H2O into one of 30 polyhedra, according to closest distance sorting
# -> for each lattice site, show number of molecules included in voronoi polyhedron 
# v01: can only count number of molecule inclusion (not total net charge per polyhedron yet)
# 4 arguments : grofile COMfile ASCIItraj output
##CAUTION : COMfile should be indexed according to sigma phase lattice site indices

import math
import sys
import numpy
import timeit

#dictionary for atomic/united atomic masses.

def gro_atomtype(str): #own needed function for splitting of .gro format.
  splitstr=str[10:15]
  splitstr=splitstr.replace(" ","")
  return splitstr

def gro_xyzsplit(str): #own needed function for splitting of .gro format. V05.
  splitstr=[str[20:28],str[28:36],str[36:44]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  posdata=numpy.array([float(x) for x in splitstr])
  return posdata

def ascii_trajsplit(str): #reading coord from traj ASCII. box size cannot be read by this.
  splitstr=[str[16:28],str[29:42],str[43:56]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  posdata=numpy.array([float(x) for x in splitstr])
  return posdata

def latticeassign(r1,box,lat): #distance calc, find closest one among lattice points
  #lat : 2D array of xyz values of all lattice points
  mindist=box[0]
  minindex=0
  for i in range(lat.shape[0]): #i: index for lattice point
    r=lat[i]-r1
    for j in range(3):
      if r[j]>(box[j]/2.0):
        r[j]-=box[j]
      elif r[j]<(-box[j]/2.0):
        r[j]+=box[j]
    dist=numpy.sqrt(numpy.dot(r,r))
    if dist<mindist:
      mindist=dist
      minindex=i
  minindex+=1 #to make the index starts from 1(not from zero)
  return minindex

def main():
  start_time=timeit.default_timer()
  #Load files
  grofile = open(sys.argv[1],'r')
  comfile = open(sys.argv[2],'r')
  trjfile = open(sys.argv[3],'r')
  outfile = open(sys.argv[4],'w') #output file for number of molecule inclusion

  #load grofile -> get atomtype-atomindex relationship first
  entry=[]
  lindex=0
  for line in grofile:
    if lindex !=0:
      if lindex==1:
        totnatom=int(line)
      elif lindex>=2 and lindex<=1+totnatom:
        split=gro_atomtype(line)
        entry.append(split)
    lindex+=1

  #load COMfile, ASCII traj simultaneously (big process!)
  lindex,sindex,effsindex=0,0,0
  comcrd=numpy.empty((0,3),float)
  waterc=numpy.zeros(30,float) #water inclusion count
  tmac=numpy.zeros(30,float) #TMA inclusion count
  for line in comfile:
    if lindex==0:
      split=line.split()
      grostep=int(split[3])
      print('read grostep {}'.format(grostep))
    if lindex !=0:
      if lindex==1:
        totlatp=int(line)
      elif lindex>=2 and lindex<=1+totlatp:
        split=gro_xyzsplit(line)
        comcrd=numpy.vstack((comcrd,split))
      elif lindex==2+totlatp:
        split=line.split()
        box=numpy.array([float(x) for x in split])
        #conclusion for 1 step, and initialization again
        #inside here, ASCII traj loading
        lindex2=0
        while lindex2 < 7+totnatom:
          trajline=trjfile.readline()
          if lindex2==0:
            sindex+=1
            split=trajline.split()
            split[2]=split[2].replace(":","")
            trjstep=int(split[2])
            if grostep==trjstep:
              effsindex+=1
              print('grostep&trjstep matched. sindex {} step {}'.format(sindex,grostep))
          if lindex2==3 and grostep==trjstep: #box vector x
            splitstr=[trajline[18:30]]
            for i in range(len(splitstr)):
              splitstr[i]=splitstr[i].replace(" ","")
            box=numpy.array([float(splitstr[0])])
          elif lindex2==4 and grostep==trjstep: #box vector y
            splitstr=[trajline[31:44]]
            for i in range(len(splitstr)):
              splitstr[i]=splitstr[i].replace(" ","")
            box=numpy.append(box,[float(splitstr[0])])
          elif lindex2==5 and grostep==trjstep: #box vector z
            splitstr=[trajline[45:58]]
            for i in range(len(splitstr)):
              splitstr[i]=splitstr[i].replace(" ","")
            box=numpy.append(box,[float(splitstr[0])])
          if lindex2>=7 and lindex2<=6+totnatom and grostep==trjstep:
            aindex=lindex2-7
            if entry[aindex]=='OW': #water : oxygen atom
              pos=ascii_trajsplit(trajline)
              addr=latticeassign(pos,box,comcrd)
              waterc[addr-1]+=1
            elif entry[aindex]=='N+': #TMA : nitrogen atom
              pos=ascii_trajsplit(trajline)
              addr=latticeassign(pos,box,comcrd)
              tmac[addr-1]+=1
          if lindex2==6+totnatom and grostep>trjstep:
            lindex2=-1
          lindex2+=1
        if (sindex%10)==0:
          elapsed=timeit.default_timer() - start_time
          print('step {} complete time {}'.format(sindex,elapsed))
        comcrd=numpy.empty((0,3),float)
        lindex=-1
    lindex+=1

  #printing section
  waterc=waterc/effsindex
  tmac=tmac/effsindex

  outfile.write('lattp# H2Ocount TMAcount (per snapshot)\n')
  for i in range(30):
    outfile.write('{:3} {:7.3f} {:7.3f}\n'.format(i+1,waterc[i],tmac[i]))

  grofile.close()
  comfile.close()
  trjfile.close()
  outfile.close()

if __name__ == "__main__": main()
