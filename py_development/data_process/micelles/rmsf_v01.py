#!/home/kjeong23/softwares/bin/python3.4
# micelle COM position rms fluctuation calculation
#<algorithm>
#starting from reindexed COM trajectory. reading -> should apply 'smart pbc reading'.
#consider 2 subsequent steps -> if displacement >= box/2, transpose the later step.
#-> collect along the whole COMtraj:get average position. -> 2nd reading(apply smart reading again)
#-> calculate deviation from ave -> get rmsf components, norm (sqrt (dx^2+dy^2+dz^2))
#2 arguments: input(reindexed COM traj) output(table of site, xyz components of rms position fluct, norm)

import math
import sys
import numpy

def gro_xyzsplit(str): #own needed function for splitting of .gro format. V05.
  splitstr=[str[20:28],str[28:36],str[36:44]]
  for i in range(len(splitstr)):
    splitstr[i]=splitstr[i].replace(" ","")
  for i in range(0,3):
    splitstr[i]=float(splitstr[i])
  return splitstr

def smartcrd(crd1,crd2,box,nmic): #'smart pbc reading'
  newcrd=crd2
  for i in range(nmic):
    for j in range(3):
      if crd2[i][j]-crd1[i][j] >= (box[j]/2.0):
        newcrd[i][j]-=box[j]
      elif crd2[i][j]-crd1[i][j] <= (-box[j]/2.0):
        newcrd[i][j]+=box[j]
  return newcrd

#main fxn
def main():
 
  #Load input files
  trjfile = open(sys.argv[1],'r') #reindexed com trajectory file
  outfile = open(sys.argv[2],'w') #output file for rms position fluctuation

  #start the loop of 'COM trajectory reading'(1st: getting ave)
  sindex,lindex=0,0
  crd1=numpy.empty((0,3),float)
  crd2=numpy.empty((0,3),float)
  for line in trjfile:
    if lindex!=0:
      if lindex==1:
        nmic=int(line)
      elif lindex>=2 and lindex<=1+nmic:
        split=gro_xyzsplit(line)
        if sindex==0: #initial step
          crd1=numpy.vstack((crd1,split))
        else:
          crd2=numpy.vstack((crd2,split))
      elif lindex==2+nmic:
        split=line.split()
        if sindex==0: #initial step
          box1=numpy.array([float(x) for x in split])
        else:
          box2=numpy.array([float(x) for x in split])

        if sindex==0:
          avecrd=crd1
        else:
          crd2=smartcrd(crd1,crd2,box2,nmic)
          avecrd+=crd2
          crd1=crd2
          crd2=numpy.empty((0,3),float)
        #initialization
        if (sindex%50)==0:
          print('loop 1 step {} complete'.format(sindex))
        sindex+=1
        lindex=-1
    lindex+=1
  avecrd/=sindex

  #Load input file again for 2nd loop, 2nd loop:getting STDEV
  trjfile = open(sys.argv[1],'r') 
  sindex,lindex=0,0
  crd1=numpy.empty((0,3),float)
  crd2=numpy.empty((0,3),float)
  stdev=numpy.zeros((nmic,3),float)
  norm=numpy.zeros(nmic,float)

  for line in trjfile:
    if lindex!=0:
      if lindex==1:
        nmic=int(line)
      elif lindex>=2 and lindex<=1+nmic:                                                              
        split=gro_xyzsplit(line)
        if sindex==0: #initial step
          crd1=numpy.vstack((crd1,split))                                                             
        else:
          crd2=numpy.vstack((crd2,split))                                                             
      elif lindex==2+nmic:
        split=line.split()
        if sindex==0: #initial step
          box1=numpy.array([float(x) for x in split])                                                 
        else:
          box2=numpy.array([float(x) for x in split])                                                 
      
        if sindex==0:
          onedev=crd1-avecrd #deviation
        else:                                                                                         
          crd2=smartcrd(crd1,crd2,box2,nmic)                                                          
          onedev=crd2-avecrd
          crd1=crd2                                                                                   
          crd2=numpy.empty((0,3),float)
        onedev=onedev*onedev
        stdev+=onedev #first adding (dev)^2 of 1 step 
        #initialization
        if (sindex%50)==0:                                                                            
          print('loop 2 step {} complete'.format(sindex))                                             
        sindex+=1
        lindex=-1
    lindex+=1
  norm=stdev[:,0]+stdev[:,1]+stdev[:,2]
  stdev/=(sindex-1) #variance
  norm/=(sindex-1) #variance

  stdev=numpy.sqrt(stdev)
  norm=numpy.sqrt(norm)

  #printing section
  for i in range(nmic):
    outfile.write('{:2} {:8.5f} {:8.5f} {:8.5f} {:8.5f}\n'.format(i+1,stdev[i][0],stdev[i][1],stdev[i][2],norm[i])) 

  trjfile.close()
  outfile.close()

if __name__ == "__main__": main()

