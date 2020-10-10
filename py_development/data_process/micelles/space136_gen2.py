#!/home/kjeong23/softwares/bin/python3.4
# program for space group operation for P42/mnm (no.136)
# to generate sigma phase

import math
import sys

#basis points (5 for sigam phase alloy. ref:Inorg. Chem. 2013, 52, 3674−3686)
#use another ref: Science 15 October 2010: vol. 330 no. 6002 pp. 349-353  -- CANCELLED
basis=[[0,0,0],[0.399,0.399,0],[0.464,0.131,0],[0.741,0.066,0],[0.187,0.187,0.251]]
#basis=
entry=[]
la=[13.34,13.34,7.01]

outfile = open(sys.argv[1],'w') #output file for sigma phase example coord

def space136(r,flag):#space group operation for P42/mnm (no.136)
  x,y,z=r[0],r[1],r[2]
  if flag==1:
    a=[x,y,z]
  elif flag==2:
    a=[-x,-y,z]
  elif flag==3:
    a=[-y+0.5,x+0.5,z+0.5]
  elif flag==4:
    a=[y+0.5,-x+0.5,z+0.5]
  elif flag==5:
    a=[-x+0.5,y+0.5,-z+0.5]
  elif flag==6:
    a=[x+0.5,-y+0.5,-z+0.5]
  elif flag==7:
    a=[y,x,-z]
  elif flag==8:
    a=[-y,-x,-z]
  elif flag==9:
    a=[-x,-y,-z]
  elif flag==10:
    a=[x,y,-z]
  elif flag==11:
    a=[y+0.5,-x+0.5,-z+0.5]
  elif flag==12:
    a=[-y+0.5,x+0.5,-z+0.5]
  elif flag==13:
    a=[x+0.5,-y+0.5,z+0.5]
  elif flag==14:
    a=[-x+0.5,y+0.5,z+0.5]
  elif flag==15:
    a=[-y,-x,z]
  elif flag==16:
    a=[y,x,z]
  a=[a[0]*la[0],a[1]*la[1],a[2]*la[2]]
  for i in range(3):
    if a[i]<0:
      a[i]+=la[i]
    elif a[i]>la[i]:
      a[i]-=la[i]
  a=[round(a[0],3),round(a[1],3),round(a[2],3)]
  return a

for r in basis:
  for i in range(1,17):
    rnew=space136(r,i)
    if rnew in entry:
      det=True 
    else:
      entry.append(rnew)

#output gro file printing section
outindex=0
outfile.write('Sigma phase alloy example coord \n')
outfile.write(' 30\n')
for r in entry:
  outindex+=1
  outfile.write('{:5}{:5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}\n'.format(outindex,'SIG','OM',outindex,r[0],r[1],r[2]))
outfile.write('         {}  {}  {}'.format(la[0],la[1],la[2]))


