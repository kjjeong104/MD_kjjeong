#!/home/kjeong23/softwares/bin/python3.4
# prototype program for micelle detection from coordinate file
# algorithm: get coord&info -> starts loop. -> distance test(consider pbc)
# -> include -> further test again&again -> full 1micelle list -> store
# -> start looking for another micelle
# v05 : able to read trajectory ascii file for statistical treatment
# ** ascii trajectory file: requires pbc -whole treatment
#(v04: size, asphericity). v041 : size-> Rg. aspher -> normalization check

import math
import sys
import numpy

def gro_minorsplit(str): #own needed function for splitting of .gro format. V05.
  splitstr=[str[0:10],str[10:15]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  return splitstr

def ascii_trajsplit(str): #reading coord from traj ASCII. box size cannot be read by this.
  splitstr=[str[16:28],str[29:42],str[43:56]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  for i in range(0,3):
    splitstr[i]=float(splitstr[i])
  return splitstr

#Part to load coordinate file, split and put into array. V05.
grofile = open(sys.argv[1],'r')
trjfile = open(sys.argv[2],'r')
outfile = open(sys.argv[3],'w') #output file for result summary
comfile = open(sys.argv[4],'w') #output file for reduced expression of micelle
minagg=int(input("minimum aggregation number to be recognized as a micelle?recomm:6\n"))
thr=float(input("threshold for micelle detection by tail position? recomm:0.9\n"))
surfstr='DEP'
entry,tendlist=[],[]
nsurf,lindex=0,0
for line in grofile:
  if lindex !=0:
    if lindex==1:
      totnatom=int(line)
      print(totnatom)
    elif lindex>=2 and lindex<=1+totnatom:
      split=gro_minorsplit(line)
      entry.append(split)
      if split[1]=='C1':
        nsurf+=1
        tendlist.append(lindex-2) #atomindex(from 0) list for tailends
  lindex+=1

#dictionary for atomic/united atomic masses. V05.
amass={'C1':15.035, 'C2':14.027, 'C3': 14.027, 'C4': 14.027, 'C5': 14.027, \
'C6': 14.027, 'C7': 14.027, 'C8': 14.027, 'C9': 14.027, 'CA': 14.027, \
'P1': 30.9738, 'OD1': 15.9994, 'OD2': 15.9994, 'OD3': 15.9994}

def distdet(r1,r2,thr,box): #distance checking function. returns true for detection
  #v04 upgrade: when thr is set to zero, then serve other function(returns dist itself)
  r=r2-r1
  for i in range(3):
    if r[i]>(box[i]/2.0):
      r[i]-=box[i]
    elif r[i]<(-box[i]/2.0):
      r[i]+=box[i]
  dist=numpy.sqrt(numpy.dot(r,r))
  if thr==0:
    return dist
  else:
    if dist<=thr:
      det=True 
    else:
      det=False
    return det

def lookup1d(n,array):#1d array lookup for a particular number
  det=False
  if n in array:
    det=True
  return det

def lookup2d(n,array):#2d array lookup for a particular number
  det=False
  for row in array:
    if n in row:
      det=True
  return det

def morphol(list,box,tempcrd): 
  #fxn to calculate morphology(COM,diameter,asphericity) of an atom group, considering pbc.
  #differently from v02, totcrd of all atoms are considered!
  #masses should be also considered by using dictionary!
  c=numpy.array([0.0, 0.0, 0.0])
  totmass=0.0
  #refinement of coordinate for wall-crossing micelles
  tends,moveends,oxys=[],[],[]
   #array for tailend atom indexes, 'marked molecules' for moving (-box), oxygens
  for i in list:
    if entry[i][1]=='C1':
      tends.append(i) #record tail ends. wall-crossing examination will be done
    elif entry[i][1]=='OD1' or entry[i][1]=='OD2' or entry[i][1]=='OD3':
      oxys.append(i) #record oxygen atoms. Micelle diameter and asphericity will be calcd.
  for i in tends: #wall-crossing test: done for only tail ends.
    for j in tends:
      for k in range(3):
        if (tempcrd[j][k]-tempcrd[i][k])>(box[k]/2.0): #if j is cut by (1b)edge
          if [entry[j][0],k] not in moveends:
            moveends.append([entry[j][0],k])
        elif (tempcrd[i][k]-tempcrd[j][k])>(box[k]/2.0): #if i is cut by (1b)edge
          if [entry[i][0],k] not in moveends:
            moveends.append([entry[i][0],k])

  for x in moveends:
    for y in list:
      if entry[y][0]==x[0]:
        tempcrd[y][x[1]]-=box[x[1]] #negative shift for wall-crossing micelle atoms
##part 1 : calculating COM
  for x in list:
    m1=amass[entry[x][1]] #corresponding atom's mass
    c+=m1*tempcrd[x]
    totmass+=m1
  c/=totmass
##part 2: micelle radius of gyration.(using only oxygens.)
##part 3: asphericity (using only oxygens. So, calculate oxygens' COM.)
 #definition T_mn = 1/2N^2 sum(i)sum(j)(r(i)_m-r(j)_m)(r(i)_n-r(j)_n)
 #or T_mn = 1/N sum(i) r(i)_m * r(i)_n when origin is set to COM
  oc=numpy.array([0.0, 0.0, 0.0])
  for x in oxys:
    oc+=tempcrd[x]
  oc/=len(oxys)
  T=numpy.zeros((3,3)) #gyration tensor. need to be 3*3 array
  for m in range(3): # dimension index 1
    for n in range(3): # dimension index 2
      for x in oxys:
        T[m,n]+=(tempcrd[x][m]-oc[m])*(tempcrd[x][n]-oc[n])
      T[m,n]/=len(oxys)
  #diagonalize T (Also, Rg^2 = lmda_x^2 + lmda_y^2 + lmda_z^2)
  e_values, e_vectors = numpy.linalg.eig(T)
  e_values=numpy.sort(e_values)
  rg2=e_values[0]+e_values[1]+e_values[2]
  rg=numpy.sqrt(rg2)
  e10,e20,e21 = e_values[1]-e_values[0],e_values[2]-e_values[0],e_values[2]-e_values[1]
  asp=(e10*e10+e20*e20+e21*e21)/(2.0*rg2*rg2) #asphericity

 #printing/returning section
  rg=round(rg,4)
  c=numpy.array([round(c[0],3),round(c[1],3),round(c[2],3)])
  d=numpy.array([c[0],c[1],c[2],rg,asp])
  return d

#main fxn:main administrator of algorithm progression
def main():
  #Whole Loop of coordinate reading!
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
    lindex+=1

    if lindex>=7+totnatom: #conclusion for 1 step, and initialization(a big process!)
      totmic,already=[],[]
      #print('step {}'.format(sindex)) #debug
      #print(totcrd) #debug
      for i in range(nsurf):#If not contained: starting a new micelle
        onemic=[]
        if lookup2d(i,totmic)==False:
          onemic.append(i)
          search=1
        while search==1:#starting of repetitive searching algorithm
          search=0
          for k in onemic:#stops when increase of number in onemic stopped
            if lookup1d(k,already)==False:
              for j in range(nsurf):
                if j!=k and lookup2d(j,totmic)==False and lookup1d(j,onemic)==False:
                  if(distdet(totcrd[tendlist[k]],totcrd[tendlist[j]],thr,box)): #detection
                    onemic.append(j)
                    search=1
              already.append(k)
        if len(onemic)!=0:
          onemic.sort()
          totmic.append(onemic)
        already=[]
     
      #list refinement section
      refmic,hom=[],[]
      unagg=0
      for row in totmic:
        if len(row)>=minagg:
          refmic.append(row)
        else:
          hom.append(row)
          unagg+=len(row)

      #com reduction, morphology calculation 
      coms,dias,asps=numpy.empty((0,3),float),numpy.array([]),numpy.array([])
      for row in refmic: #tailend member list of 1 micelle
        onemic=[]
        for i in row:
          for j in range(totnatom):
            if ((str(i+1)+surfstr)==entry[j][0]): #search for specific molID
              onemic.append(j) #append atomindex
        d=morphol(onemic,box,totcrd) #total data array
        c=numpy.array([d[0],d[1],d[2]])
        coms=numpy.vstack((coms,c))
        dias=numpy.append(dias,d[3])
        asps=numpy.append(asps,d[4])

      #printing section
      micindex=0
      outfile.write('step# {} numofmcls {}\n'.format(sindex,len(refmic)))
      for row in refmic:
        micindex+=1
        outfile.write('mic# {} N= {} r= {} A= {} members: {}\n'.format(micindex,len(row),dias[micindex-1],round(asps[micindex-1],3),row))
      outfile.write('{} Unaggregated molecules: {}\n'.format(unagg,hom))

      #com output writing section
      outindex=0
      comfile.write('Centers of micelles: step {}\n'.format(sindex))
      comfile.write(' {}\n'.format(micindex))
      for row in coms:
        outindex+=1
        comfile.write('{:5}{:5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}\n'.format(outindex,'COM','X',outindex,row[0],row[1],row[2]))
      comfile.write('      {}  {}  {}\n'.format(box[0],box[1],box[2]))

      #initialization statements(at last)
      box=numpy.array([])
      totcrd=numpy.empty((0,3),float)
      lindex=0
      sindex+=1 
  
  grofile.close()
  trjfile.close()
  outfile.close()
  comfile.close()

if __name__ == "__main__": main()

