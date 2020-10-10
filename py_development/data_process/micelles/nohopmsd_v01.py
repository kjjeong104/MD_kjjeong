#!/home/kjeong23/softwares/bin/python3.4
# program to calculate hoppingless MSD of surfactants
# algorithm: track 1 surfactant. Combine coords info with micelle inclusion info.
# -> record r0 and COM. -> calc square displacement until it hops.(consider PBC)
# in every step of COM SD calc, pile up in the bins, adding up "counts"(w.r.t. t)
# to compare statistical weight of all dr^2(t). -> calc MSD(t).
# -> display MSD(t), count(t).
# ** ascii trajectory file: requires pbc -whole treatment!!
#inputs: grofile, traj ASCII file, morfile(lat assigned)(to get site info) 
#output: file of {t, nonhop_MSD(t), count(t)}

import math
import sys
import numpy
import timeit

#atomic mass dictionary
amass={'C1':15.035, 'C2':14.027, 'C3': 14.027, 'C4': 14.027, 'C5': 14.027, \
'C6': 14.027, 'C7': 14.027, 'C8': 14.027, 'C9': 14.027, 'CA': 14.027, \
'P1': 30.9738, 'OD1': 15.9994, 'OD2': 15.9994, 'OD3': 15.9994, \
'N+': 14.0067, 'Me1': 15.035, 'Me2': 15.035, 'Me3': 15.035, 'Me4': 15.035, \
'OW': 15.9994, 'HW1': 1.008, 'HW2': 1.008}

def gro_minorsplit(str): #own needed function for splitting of .gro format. V05.
  splitstr=[str[0:10],str[10:15]] #molecule index&name, atomname
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  return splitstr

def ascii_trajsplit(str): #reading coord from traj ASCII. box size cannot be read by this.
  splitstr=[str[16:28],str[29:42],str[43:56]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  posdata=numpy.array([float(x) for x in splitstr])
  return posdata

def comcalc(alist,crd): #calculation of COM of a molecule
  compos=numpy.zeros(3)
  totmass,i=0.0,0
  for x in alist: #alist : list of atomnames. #crd: n*3 array of coords
    m1=amass[x]
    totmass+=m1
    compos+=m1*crd[i]
    i+=1
  compos/=totmass
  return compos

def pbcsd(ri,rf,box,crcount): #calc of square displacement with considering pbc
  #crcount : 1*3 array of box wall crossing count of that specific molecule
  rf+=box*crcount
  r=rf-ri
  sd=numpy.dot(r,r)
  return sd

#main fxn
def main():
  #Part to load coordinate file, morphology info file
  grofile = open(sys.argv[1],'r')
  trjfile = open(sys.argv[2],'r')
  morfile = open(sys.argv[3],'r')
  outfile = open(sys.argv[4],'w')

  initt=float(input("Initial time of this trajectory (in ns)?\n"))
  fint=float(input("Final time of the whole trajectory in ns? ex) 499.96 \n"))
  ssize=float(input("Timestep between snapshots (in ns)? ex) 0.04\n"))
  surfstr='DEP'
  nsurf,lindex,entry=0,0,[]
  start_time=timeit.default_timer()

  #reading grofile to match molecule,atom index
  for line in grofile:
    if lindex !=0:
      if lindex==1:
        totnatom=int(line)
      elif lindex>=2 and lindex<=1+totnatom:
        split=gro_minorsplit(line)
        entry.append(split)
        if split[1]=='C1':
          nsurf+=1
    lindex+=1
  print(nsurf) #debug

  #prepare bins
  nstep=int((fint-initt)/ssize)
  sdbin,count=numpy.zeros(nstep),numpy.zeros(nstep) #bin & stat weight for dr^2(t) ave
  #rp : previous step position of COM
  r0,rp,rt=numpy.zeros((nsurf,3)),numpy.zeros((nsurf,3)),numpy.zeros((nsurf,3))
  cross=numpy.zeros((nsurf,3)) #count for boundary crossing of each surfactant

  #loop of reading micelle inclusion data & coord steps (including decision of r0)
  lindex,sindex=0,0
  totcrd=numpy.empty((0,3),float)
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
    elif lindex>=7 and lindex<=6+totnatom:
      split=ascii_trajsplit(line)
      totcrd=numpy.vstack((totcrd,split))
    if lindex==6+totnatom: #conclusion for 1 step, initialize
      tempcrd=numpy.empty((0,3),float) #coord set for 1molecule
      alist=[]
      aindex=0
      #calc of COMs
      for row in entry: #search molecule index & atom names
        molnum=int(row[0].replace(surfstr,""))
        alist.append(row[1])
        tempcrd=numpy.vstack((tempcrd,totcrd[aindex]))
        if row[1]=='OD3': #final atom of a molecule
          #print('{} {}'.format(row[0],tempcrd))
          if sindex==0:
             r0[molnum-1]=comcalc(alist,tempcrd)
             rp[molnum-1]=r0[molnum-1]
          else:
             rt[molnum-1]=comcalc(alist,tempcrd)
             #crossing test
             for k in range(3):
               if rt[molnum-1][k]-rp[molnum-1][k]>(box[k]/2.0): #leftside.crossing
                 cross[molnum-1][k]-=1 #should transpose final position (-)
               elif rt[molnum-1][k]-rp[molnum-1][k]<(-box[k]/2.0): #rightside.crossing
                 cross[molnum-1][k]+=1 #should transpose final position (+)
             rp[molnum-1]=rt[molnum-1] #after crossing test, store position
          #initialize single-molecule variables
          alist=[]
          tempcrd=numpy.empty((0,3),float)
        aindex+=1
      #after calc of all COMs
      for j in range(nsurf):
        if sindex>0:
          #print(r0[j],rt[j])
          sd=pbcsd(r0[j],rt[j],box,cross[j])
          sdbin[sindex-1]+=sd
          count[sindex-1]+=1
      print(cross)
      totcrd=numpy.empty((0,3),float)
      sindex+=1
      lindex=-1
      elapsed=timeit.default_timer() - start_time
      print('going to step {} time {:8.4f}'.format(sindex,elapsed))
    lindex+=1

  #printing section
  sdbin/=count
  outfile.write('{:8.4f} {:8.4f} {:8.4f}\n'.format(0,0,0))
  for i in range(1,nstep+1):
    outfile.write('{:8.4f} {:8.4f} {:8.4f}\n'.format(ssize*i,sdbin[i-1],count[i-1]))

  grofile.close()
  trjfile.close()
  morfile.close()
  outfile.close()

if __name__ == "__main__": main()

