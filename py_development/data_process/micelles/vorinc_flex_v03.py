#!/home/kjeong23/softwares/bin/python3.4
# Voronoi polyhedron counterion/water inclusion counting program (for any micellar system)
# algorithm: Don't actually calculate Voronoi polyhedron! 
# get atomindex info from grofile -> read COMfile(1step) and ASCII traj(1step) simultaneously
# -> classify each MeN, H2O into one of (detected#) polyhedra, according to closest distance sorting
# -> for each micelle(tracked by other script), show #of molecules included in voronoi polyhedron 
# v02: can calculate total mass and total charge in each voronoi cell (not the density yet)
# v021: can determine the number of lattice sites depending on the input
# v03: now for every snapshot, integrates morphology info and voronoi polyhedron inclusion
# doesn't get user input for number of lattice sites. It reads from COMfile
# 5 arguments : grofile ASCIItraj COMfile morfile output
# caution : for COMfile and morfile, use raw files without reindexing.
# this is tailored for non-characterized micellar systems.

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

def txt_morsplit(str): #reading morphology summary file.
  splitstr=[str[5:8],str[12:15],str[18:27],str[30:39],str[47:56],str[56:65],str[65:74]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  mordata=numpy.array([float(x) for x in splitstr])
  return mordata

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
  trjfile = open(sys.argv[2],'r')
  comfile = open(sys.argv[3],'r')
  morfile = open(sys.argv[4],'r')
  outfile = open(sys.argv[5],'w') #output file for number of molecule inclusion

  #totlatp=int(input("Designated # of sites for Voronoi analysis? ex) 33 \n"))
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

  #load COMfile, ASCII traj, morfile simultaneously (big process!)
  lindex,sindex=0,0
  comcrd=numpy.empty((0,3),float)
#  waterc,tmac,massc,chgc=numpy.zeros(totlatp,float),numpy.zeros(totlatp,float),numpy.zeros(totlatp,float),numpy.zeros(totlatp,float) 
  for line in comfile:
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
        #inside here, ASCII traj loading and morfile loading
        lindex2=0
        waterc,tmac,massc,chgc=numpy.zeros(totlatp,int),numpy.zeros(totlatp,int),numpy.zeros(totlatp,float),numpy.zeros(totlatp,float)
        while lindex2 < 7+totnatom:
          trajline=trjfile.readline()
          if lindex2==3: #box vector x
            splitstr=[trajline[18:30]]
            for i in range(len(splitstr)):
              splitstr[i]=splitstr[i].replace(" ","")
            box=numpy.array([float(splitstr[0])])
          elif lindex2==4: #box vector y
            splitstr=[trajline[31:44]]
            for i in range(len(splitstr)):
              splitstr[i]=splitstr[i].replace(" ","")
            box=numpy.append(box,[float(splitstr[0])])
          elif lindex2==5: #box vector z
            splitstr=[trajline[45:58]]
            for i in range(len(splitstr)):
              splitstr[i]=splitstr[i].replace(" ","")
            box=numpy.append(box,[float(splitstr[0])])
          if lindex2>=7 and lindex2<=6+totnatom:
            aindex=lindex2-7
            pos=ascii_trajsplit(trajline)
            addr=latticeassign(pos,box,comcrd)
            massc[addr-1]+=amass[entry[aindex]]
            chgc[addr-1]+=acharge[entry[aindex]]
            if entry[aindex]=='OW': #water : oxygen atom
              waterc[addr-1]+=1
            elif entry[aindex]=='N+': #TMA : nitrogen atom
              tmac[addr-1]+=1
          lindex2+=1
        #morfile loading and output printing section
        lindex3=0
        while lindex3 < 2+totlatp:
          morline=morfile.readline()
          if lindex3==0:
            split=morline.split()
            outfile.write('mic# agN    Rg       Asp     sv2a     sv2b     sv2c \
  waterc  TMAc totmass totchg   step {} numofmcls {}\n'.format(sindex,split[3]))
          elif lindex3>=1 and lindex3 <= totlatp:
            split=txt_morsplit(morline)
            i=lindex3-1
            outfile.write('{:3} {:3} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f} {:5} {:5} {:8.4f} {:8.4f}\n'\
.format(int(split[0]),int(split[1]),split[2],split[3],split[4],split[5],split[6],waterc[i],tmac[i],massc[i],chgc[i]))
          lindex3+=1
        if (sindex%10)==0:
          elapsed=timeit.default_timer() - start_time
          print('step {} complete time {}'.format(sindex,elapsed))
        comcrd=numpy.empty((0,3),float)
        sindex+=1
        lindex=-1
    lindex+=1
  #printing section
#  waterc=waterc/effsindex
#  tmac=tmac/effsindex
#  massc=massc/effsindex
#  chgc=chgc/effsindex

#  outfile.write('lattp# H2Ocount TMAcount totmass totchg (per snapshot)\n')
#  for i in range(8):
#    outfile.write('{:3} {:7.3f} {:7.3f} {:8.4f} {:8.4f}\n'.format(i+1,waterc[i],tmac[i],massc[i],chgc[i]))
  grofile.close()
  comfile.close()
  trjfile.close()
  morfile.close()
  outfile.close()

if __name__ == "__main__": main()
