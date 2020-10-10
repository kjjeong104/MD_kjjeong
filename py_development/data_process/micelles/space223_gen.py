#!/home/kjeong23/softwares/bin/python3.4
# program for space group operation for Pm 3bar n (no.223)
# to generate A15 morphology

import math
import sys

#basis points (2 for A15 phase alloy. 
#ref: Graef, Marc De; McHenry, Michael E (2007). Structure of materials: an introduction to crystallography, diffraction and symmetry. pp. 518–521.)
basis=[[0,0,0],[0.250,0,0.500]]
entry=[]
la=[6.935,6.935,6.935]

outfile = open(sys.argv[1],'w') #output file for sigma phase example coord

def space223(r,flag):#space group operation for Pm 3bar n (no.223)
  x,y,z=r[0],r[1],r[2]
  if flag==1:
    a=[x,y,z]
  elif flag==2:
    a=[x,-y,-z]
  elif flag==3:
    a=[-x,y,-z]
  elif flag==4:
    a=[-x,-y,z]
  elif flag==5:
    a=[z,x,y]
  elif flag==6:
    a=[-z,-x,y]
  elif flag==7:
    a=[z,-x,-y]
  elif flag==8:
    a=[-z,x,-y]
  elif flag==9:
    a=[y,z,x]
  elif flag==10:
    a=[-y,z,-x]
  elif flag==11:
    a=[-y,-z,x]
  elif flag==12:
    a=[y,-z,-x]
  elif flag==13:
    a=[x+0.5,-z+0.5,y+0.5]
  elif flag==14:
    a=[x+0.5,z+0.5,-y+0.5]
  elif flag==15:
    a=[-x+0.5,-z+0.5,-y+0.5]
  elif flag==16:
    a=[-x+0.5,z+0.5,y+0.5]
  elif flag==17:
    a=[0.5+z,0.5+y,0.5-x]
  elif flag==18:
    a=[0.5-z,0.5+y,0.5+x]
  elif flag==19:
    a=[0.5-z,0.5-y,0.5-x]
  elif flag==20:
    a=[0.5+z,0.5-y,0.5+x]
  elif flag==21:
    a=[0.5-y,0.5+x,0.5+z]
  elif flag==22:
    a=[0.5+y,0.5-x,0.5+z]
  elif flag==23:
    a=[0.5-y,0.5-x,0.5-z]
  elif flag==24:
    a=[0.5+y,0.5+x,0.5-z]
  elif flag==25:
    a=[-x,-y,-z]
  elif flag==26:
    a=[-x,y,z]
  elif flag==27:
    a=[x,-y,z]
  elif flag==28:
    a=[x,y,-z]
  elif flag==29:
    a=[-z,-x,-y]
  elif flag==30:
    a=[z,x,-y]
  elif flag==31:
    a=[-z,x,y]
  elif flag==32:
    a=[z,-x,y]
  elif flag==33:
    a=[-y,-z,-x]
  elif flag==34:
    a=[y,-z,x]
  elif flag==35:
    a=[y,z,-x]
  elif flag==36:
    a=[-y,z,x]
  elif flag==37:
    a=[0.5-x,0.5+z,0.5-y]
  elif flag==38:
    a=[0.5-x,0.5-z,0.5+y]
  elif flag==39:
    a=[0.5+x,0.5+z,0.5+y]
  elif flag==40:
    a=[0.5+x,0.5-z,0.5-y]
  elif flag==41:
    a=[0.5-z,0.5-y,0.5+x]
  elif flag==42:
    a=[0.5+z,0.5-y,0.5-x]
  elif flag==43:
    a=[0.5+z,0.5+y,0.5+x]
  elif flag==44:
    a=[0.5-z,0.5+y,0.5-x]
  elif flag==45:
    a=[0.5+y,0.5-x,0.5-z]
  elif flag==46:
    a=[0.5-y,0.5+x,0.5-z]
  elif flag==47:
    a=[0.5+y,0.5+x,0.5+z]
  elif flag==48:
    a=[0.5-y,0.5-x,0.5+z]
  a=[a[0]*la[0],a[1]*la[1],a[2]*la[2]]
  for i in range(3):
    if a[i]<0:
      a[i]+=la[i]
    elif a[i]>la[i]:
      a[i]-=la[i]
  a=[round(a[0],3),round(a[1],3),round(a[2],3)]
  return a

for r in basis:
  for i in range(1,49):
    rnew=space223(r,i)
    if rnew in entry:
      det=True 
    else:
      entry.append(rnew)

#output gro file printing section
outindex=0
outfile.write('A15 alloy example coord \n')
outfile.write(' 8\n')
for r in entry:
  outindex+=1
  outfile.write('{:5}{:5}{:5}{:5}{:8.3f}{:8.3f}{:8.3f}\n'.format(outindex,'AFI','OM',outindex,r[0],r[1],r[2]))
outfile.write('         {}  {}  {}'.format(la[0],la[1],la[2]))


