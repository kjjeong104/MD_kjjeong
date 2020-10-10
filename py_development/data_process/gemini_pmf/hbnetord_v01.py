#!/home/kjeong23/softwares/bin/python3.4
# Hydrogen bond network orderedness counter (for any water-included system)
# algorithm: define solute, then define 'layers' of H-bond sequence from solute group
# finally get statistics of #Hbond vs layer relationship
# get atomindex list of solute acceptor -> read ndxfile(1step)(3column:donorO-H-acceptor) 
# -> starting from solute acceptor, assign all water molecules' Hbond network layers
# -> count #Hbond for all water molecules, and average by layers -> output
# v01 : 1step, v02 : averaging multiple snapshots
# 3 arguments : accfile ndxfile output
# caution : the .ndx file should be obtained with '-nomerge' option in g_hbond module

import math
import sys
import numpy
import timeit

def main():
  #Load files
  accfile = open(sys.argv[1],'r')
  ndxfile = open(sys.argv[2],'r')
  outfile = open(sys.argv[3],'w') #output file for ave-number of Hbonds wrt layer

  start_time=timeit.default_timer()
  #load accfile(artificially made index file) -> get list of solute,water indices
  solute,water=[],[]
  delim_solute,delim_OW,delim_hbonds='[ solute_System ]','[ OW_System ]','[ hbonds_System ]'
  maxnl=10 #max number of layers we inspect
  flag,nl=0,1
  for line in accfile:
    if delim_solute in line:
      flag=1
    elif delim_OW in line:
      flag=2
    if flag==1:
      split=line.split()
      if split[0]!='[':
        for x in split:
          solute.append(x)        
    elif flag==2:  
      split=line.split()
      if split[0]!='[':
        for x in split:
          water.append(x)

  layers,sets,onelayer=[],[],[]
  usedw=[] #list of used water molecules for hbond
  #load ndxfile to memorize information
  for line in ndxfile:
    if delim_hbonds in line:
      flag=3
    if flag==3:
      split=line.split() #split the 3-column line into 3 entries
      if split[0]!='[':
        sets.append([split[0],split[2]]) #put pairs: donorO, acceptor atom(solute can be only acc)

  elapsed=timeit.default_timer() - start_time
  print('information loading complete time {}'.format(elapsed))

  #layering section
  layers.append(solute) #Register 0th layer(solute). layers should be 2d list. each row is a layer
  recent=solute #most recent layer
  while nl<=maxnl: #recursive loop for water layer registration
    for row in sets:
      if row[1] in recent: #if the recent layer element is acceptor
        if row[0] not in usedw:
          onelayer.append(row[0])
          usedw.append(row[0])

    layers.append(onelayer)  #register layer
    nl+=1
    recent=onelayer
    onelayer=[]

  elapsed=timeit.default_timer() - start_time
  print('layering complete time {}'.format(elapsed))

  #Hbond counting section (actually we can combine this process with the layering section)
  nhb=numpy.zeros(maxnl+1,int) #include 0th layer
  for row in sets: #count all layers in one cycle of sweeping the ndx file
    nl=0
    for onelayer in layers:
      for x in onelayer:
        if x in row:
          nhb[nl]+=1
      nl+=1

  elapsed=timeit.default_timer() - start_time
  print('Hbond counting complete time {}'.format(elapsed))

 #printing section
  outfile.write('layer# mem# tothbn avehbn\n')
  for i in range(maxnl+1):
    nmem=len(layers[i])
    outfile.write('l {:3} {:3} {:5} {:8.4f}\n'.format( i, nmem, nhb[i], nhb[i]/float(nmem) ))
  accfile.close()
  ndxfile.close()
  outfile.close()

if __name__ == "__main__": main()
