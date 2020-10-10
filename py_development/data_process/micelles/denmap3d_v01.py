#!/home/kjeong23/softwares/bin/python3.4
# 3d density map construction program (gnuplot style output : paragraphs of 2d plots)
# algorithm: put basic settings(have to determine principal axis to slice into 2d plots,
# # of bins in each dimension(for sigma phase: use multiple of 4 for z-axis)),->
# get coord&info from ASCII traj -> smart binning by dividing coord value into binsize
# -> print in 4-column format : x(binmin) y(binmin) z(binmin) ave_density(kg/m^3)
# v01: can only show average density over whole trajectory for each bin
# 3 arguments : grofile ASCIItraj output

import math
import sys
import numpy
import timeit

#dictionary for atomic/united atomic masses.
amass={'C1':15.035, 'C2':14.027, 'C3': 14.027, 'C4': 14.027, 'C5': 14.027, \
'C6': 14.027, 'C7': 14.027, 'C8': 14.027, 'C9': 14.027, 'CA': 14.027, \
'P1': 30.9738, 'OD1': 15.9994, 'OD2': 15.9994, 'OD3': 15.9994}

def gro_atomtype(str): #own needed function for splitting of .gro format.
  splitstr=str[10:15]
  splitstr=splitstr.replace(" ","")
  return splitstr

def ascii_trajsplit(str): #reading coord from traj ASCII. box size cannot be read by this.
  splitstr=[str[16:28],str[29:42],str[43:56]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  posdata=numpy.array([float(x) for x in splitstr])
  return posdata

def smartbin(box,nbin,pos): #smart binning function with considering PBC
  #examine relevance of 'pos' first (consider PBC)
  for i in range(3):
    while pos[i]>=box[i]:
      pos[i]-=box[i]
    while pos[i]<0:
      pos[i]+=box[i]
  #determine bin position
  mesh,address=numpy.zeros(3),numpy.zeros(3)
  mesh=box/nbin #binsize of x,y,z axes
  address=pos/mesh  #address: index starts from 0
  address=numpy.fix(address) #discard decimal numbers. 
  return address

def main():
  start_time=timeit.default_timer()
  #Load files
  grofile = open(sys.argv[1],'r')
  trjfile = open(sys.argv[2],'r')
  denfile = open(sys.argv[3],'w') #output file for gnuplot style density map(2d slabs)

  #settings
  prin=input("Principal axis to slice 3d density map into 2d maps? ex) z\n")
  nbinstr=input("#bins in x,y,z axes respectively? ex) 100 100 40\n")
  split=nbinstr.split()
  nbin=numpy.array([int(x) for x in split])

  #mbin3d=numpy.zeros((nbin[0],nbin[1],nbin[2]),float) #total mass bin by 3dimension(temp. for 1step)
  dbin3d=numpy.zeros((nbin[0],nbin[1],nbin[2]),float) #mass density bin by 3dimension

  #just in case (deal NPT as NVT)
  absboxstr=input("Would you like to manually set cell parameter? ex) no or 13.5 13.5 7.1\n")
  if absboxstr!='no':
    split=absboxstr.split()
    absbox=numpy.array([float(x) for x in split])

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

  #load ASCII traj(whole loop of coordinate reading)
  lindex,sindex=0,0
  for line in trjfile:
    if lindex==3: #box vector x
      splitstr=[line[18:30]]
      for i in range(len(splitstr)):
        splitstr[i]=splitstr[i].replace(" ","")
      box=numpy.array([float(splitstr[0])])
    elif lindex==4: #box vector y
      splitstr=[line[31:44]]
      for i in range(len(splitstr)):
        splitstr[i]=splitstr[i].replace(" ","")
      box=numpy.append(box,[float(splitstr[0])])
    elif lindex==5: #box vector z
      splitstr=[line[45:58]]
      for i in range(len(splitstr)):
        splitstr[i]=splitstr[i].replace(" ","")
      box=numpy.append(box,[float(splitstr[0])])
      boxsave=numpy.array([box[0],box[1],box[2]])
    elif lindex>=7 and lindex<=6+totnatom:
      pos=ascii_trajsplit(line)
      addr=smartbin(box,nbin,pos) #bin address
      aindex=lindex-7
      m1=amass[entry[aindex]] #corresponding atom's mass
      dbin3d[addr[0]][addr[1]][addr[2]]+=m1
    lindex+=1

    if lindex>=7+totnatom: #conclusion for 1 step, and initialization
      box=numpy.array([])
      lindex=0
      if sindex%50 == 0:
        elapsed=timeit.default_timer() - start_time
        print('finished step {} time {}'.format(sindex,elapsed))
      sindex+=1

  #printing section
  mesh=numpy.zeros(3)
  if absboxstr!='no':
    mesh=absbox/nbin
  else:
    mesh=boxsave/nbin
  meshvol=mesh[0]*mesh[1]*mesh[2] #volume element in nm^3
  #conversion: total cumulative amu -> kg/m^3 per snapshot
  #1 amu = 1.660539 * 10^-27 kg. 1nm^-3 = 10^27 kg^-3
  dbin3d=dbin3d*(1.660539)/(meshvol*sindex)

  if prin=='x':
    for i in range(nbin[0]):
      for j in range(nbin[1]):
        for k in range(nbin[2]):
          denfile.write('{:8.5f} {:8.5f} {:8.5f} {:8.5f}\n'.format(mesh[0]*i,mesh[1]*j,mesh[2]*k,dbin3d[i][j][k]))
        denfile.write(' \n')
      denfile.write('-----slab: x= {:8.5f}'.format(mesh[0]*i))
      denfile.write(' \n')

  elif prin=='y':
    for j in range(nbin[1]):
      for i in range(nbin[0]):
        for k in range(nbin[2]):
          denfile.write('{:8.5f} {:8.5f} {:8.5f} {:8.5f}\n'.format(mesh[0]*i,mesh[1]*j,mesh[2]*k,dbin3d[i][j][k]))
        denfile.write(' \n')
      denfile.write('-----slab: y= {:8.5f}'.format(mesh[1]*j))
      denfile.write(' \n')

  elif prin=='z':
    for k in range(nbin[2]):
      for i in range(nbin[0]):
        for j in range(nbin[1]):
          denfile.write('{:8.5f} {:8.5f} {:8.5f} {:8.5f}\n'.format(mesh[0]*i,mesh[1]*j,mesh[2]*k,dbin3d[i][j][k]))
        denfile.write(' \n')
      denfile.write('-----slab: z= {:8.5f}'.format(mesh[2]*k))
      denfile.write(' \n')

  grofile.close()
  trjfile.close()
  denfile.close()

if __name__ == "__main__": main()
