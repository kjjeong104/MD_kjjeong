#!/home/kjeong23/softwares/bin/python3.4
# Voronoi polyhedron counterion/water inclusion counting program (for any micellar system)
# algorithm: Don't actually calculate Voronoi polyhedron! 
# get atomindex info from grofile -> read COMfile(1step) and ASCII traj(1step) simultaneously
# -> classify each MeN, H2O into one of (detected#) polyhedra, according to closest distance sorting
# -> for each micelle(tracked by other script), show #of molecules included in voronoi polyhedron 
# v02: can calculate total mass and total charge in each voronoi cell (not the density yet)
# v021: can determine the number of lattice sites depending on the input
# 4 arguments : grofile COMfile ASCIItraj output
##CAUTION : COMfile should be indexed according to micelle indices with consistent tracking

import math
import sys
import numpy
import timeit

#dictionary for atomic/united atomic masses, charges.
amass={'C1':15.035, 'C2':14.027, 'C3': 14.027, 'C4': 14.027, 'C5': 14.027, \
'C6': 14.027, 'C7': 14.027, 'C8': 14.027, 'C9': 14.027, 'CA': 14.027, \
'P1': 30.9738, 'OD1': 15.9994, 'OD2': 15.9994, 'OD3': 15.9994, \
'N+': 14.0067, 'Me1': 15.035, 'Me2': 15.035, 'Me3': 15.035, 'Me4': 15.035, \
'OW': 15.9994, 'HW1': 1.008, 'HW2': 1.008}
acharge={'C1':0.0, 'C2':0.0, 'C3': 0.0, 'C4': 0.0, 'C5': 0.0, \
'C6': 0.0, 'C7': 0.0, 'C8': 0.0, 'C9': 0.0, 'CA': -0.691, \
'P1': 1.889, 'OD1': -1.066, 'OD2': -1.066, 'OD3': -1.066, \
'N+': 0.0, 'Me1': 0.25, 'Me2': 0.25, 'Me3': 0.25, 'Me4': 0.25, \
'OW': -0.82, 'HW1': 0.41, 'HW2': 0.41}

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
  #Load files
  grofile = open(sys.argv[1],'r')
  comfile = open(sys.argv[2],'r')
  trjfile = open(sys.argv[3],'r')
  outfile = open(sys.argv[4],'w') #output file for number of molecule inclusion

  totlatp=int(input("Designated # of sites for Voronoi analysis? ex) 33 \n"))
  start_time=timeit.default_timer()
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
  waterc,tmac,massc,chgc=numpy.zeros(totlatp,float),numpy.zeros(totlatp,float),numpy.zeros(totlatp,float),numpy.zeros(totlatp,float) 
  for line in comfile:
    if lindex==0:
      split=line.split()
      grostep=int(split[3])
      print('read grostep {}'.format(grostep))
    if lindex !=0:
      if lindex==1:
        junk=int(line)
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
            pos=ascii_trajsplit(trajline)
            addr=latticeassign(pos,box,comcrd)
            massc[addr-1]+=amass[entry[aindex]]
            chgc[addr-1]+=acharge[entry[aindex]]
            if entry[aindex]=='OW': #water : oxygen atom
              waterc[addr-1]+=1
            elif entry[aindex]=='N+': #TMA : nitrogen atom
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
  massc=massc/effsindex
  chgc=chgc/effsindex

  outfile.write('lattp# H2Ocount TMAcount totmass totchg (per snapshot)\n')
  for i in range(8):
    outfile.write('{:3} {:7.3f} {:7.3f} {:8.4f} {:8.4f}\n'.format(i+1,waterc[i],tmac[i],massc[i],chgc[i]))

  grofile.close()
  comfile.close()
  trjfile.close()
  outfile.close()

if __name__ == "__main__": main()
