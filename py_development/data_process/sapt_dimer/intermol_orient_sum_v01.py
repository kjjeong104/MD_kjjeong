##python script for hydrogen bond intermolecular orientation printer.
#dihedral, angle, distance.
#be careful of PBC.
import sys
import math
import mdtraj as md
#from __future__ import print_function
#from simtk.openmm.app import *
#from simtk.openmm import *
#from simtk.unit import *
#from sys import stdout
#from time import gmtime, strftime
#from datetime import datetime
import numpy as np
#once Cl-Chol swapping problem has been resolved in other code, can erase that part.

def tg_classtext(dih_value): #trans-gauche classification (cover all range of -360~360 deg for just in case)
    #in 0~180deg scale: 60 is gauche, 180 is trans. so 0~120deg gauche, 120~180deg trans.
    #in 360deg scale : 0~120 gauche, 120~240 trans, 240~360 gauche.
    #-360~360 deg covering: -360~-240 gauche, -240~-120 trans, -120~0~120 gauche, 120~240 trans, 240~360 gauche.
    if dih_value>=-360 and dih_value<-240:
      text_tg='gau'
    elif dih_value>=-240 and dih_value<-120:
      text_tg='tra'
    elif dih_value>=-120 and dih_value<120:
      text_tg='gau'
    elif dih_value>=120 and dih_value<240:
      text_tg='tra'
    elif dih_value>=240 and dih_value<=360:
      text_tg='gau'
    return text_tg

#Structure : determine type of analysis.
#SAPT-pol cholOHCl, SAPT-pol cholOHuO, FPMD cholOHCl, FPMD cholOHuO
intcode=int(input("Put code. SAPT-pol cholOHCL:0,SAPT-pol cholOHuO:1,FPMD cholOHCl:2,FPMD cholOHuO:3\n"))
outstrlist=['saptpol_cholOHCl','saptpol_cholOHuO','fpmd_cholOHCl','fpmd_cholOHuO']
inputstr=['Chol_Cl_','Urea_Chol_','Chol_Cl_','Chol_Urea_']
textlist_sub_dist1=['NCldist','chHuNdist','NCldist','chHuNdist']
textlist_sub_ang1=['OHClang','OHuOang','OHClang','OHuOang']
#textlist_sub_dih2=['null','chCOuONdih','null','chCOuONdih']
#outstr,text_sub_dist1,text_sub_ang1,text_sub_dih2=\
#outstrlist[intcode],textlist_sub_dist1[intcode],textlist_sub_ang1[intcode],textlist_sub_dih2[intcode]
outstr,text_sub_dist1,text_sub_ang1=outstrlist[intcode],textlist_sub_dist1[intcode],textlist_sub_ang1[intcode]

#atom indices list for geometry calculation (in comment: pdbfile. need to -1 from that)
#NCl SAPTpol : 1 29/ chHuN SAPTpol : 33 (2,3) /
#NCl FPMD : 16 22 /chHuN FPMD : 21 (23,25) 
mother_list_dist=np.array([ [0,28],[32,1],[15,21],[20,22] ]) #
mother_list_dist_alter=np.array([ [0,28],[32,2],[15,21],[20,24] ])
#OHCl SAPTpol: 20 21 29, OHuO SAPTpol: 32 33 1/
# OHCl FPMD:20 21 22/OHuO FPMD: 20 21 24 
mother_list_angle=np.array([ [19,20,28],[31,32,0],[19,20,21],[19,20,23] ])
#NCCO : SAPT cholOHCl: 1 14 17 20 / SAPT cholOHuO: 13 26 29 32/ 
#FPMD cholOHCl:16 13 17 20 / FPMD cholOHuO 16 13 17 20
mother_list_dih1=np.array([ [0,13,16,19],[12,25,28,31],[15,12,16,19],[15,12,16,19] ])
#null / chCOuON SAPTpol: 30 33 1 2(choose 2or3) / null / chCOuON FPMD: 17 20 24 23
#mother_list_dih2=np.array([ [0,1,2,3],[29,32,0,1],[0,1,2,3],[16,19,23,22] ])
ndimers=1000
outfile= open("sum_intgeom_"+outstr+".dat",'w')

list_distc1=np.array([mother_list_dist[intcode]]) #2-dimensionalize to fit syntax
list_distc2=np.array([mother_list_dist_alter[intcode]])
list_angle=np.array([mother_list_angle[intcode]])
list_dih1=np.array([mother_list_dih1[intcode]])
#list_dih2=np.array([mother_list_dih2[intcode]])

#loop : dimers
for i in range(ndimers):
  pdbin=inputstr[intcode]+str(i)+'.pdb'
  traj=md.load(pdbin,top=pdbin)
  topology=traj.topology
  #if intcode==2:
  #be careful. if chol cl order swapped(due to residue index), need to change calculation
  #  if 'Cl' in str(topology.atom(0)): #first atom is Cl -> chol cl swapped
  #    list_distc1=np.array([mother_list_dist[4]])
  #    list_angle=np.array([mother_list_angle[4]])
  #    list_dih1=np.array([mother_list_dih1[4]])

  sub_distc1=(md.compute_distances(traj,list_distc1)).flatten()[0]
  sub_distc2=10000000
  sub_ang1=(md.compute_angles(traj,list_angle)).flatten()[0]
  if intcode==1 or intcode==3:
    sub_distc2=(md.compute_distances(traj,list_distc2)).flatten()[0]
  #if intcode==0:
  #elif intcode==1:
  #chHuN distance: choose closer one
  #elif intcode==2: 
  #be careful. if chol cl order swapped(due to residue index), need to change calculation
  #elif intcode==3:
  #chHuN distance: choose closer one
  sub_dist1=np.minimum(sub_distc1,sub_distc2)
  main_dih=(md.compute_dihedrals(traj,list_dih1)).flatten()[0]
  #sub_dih2=(md.compute_dihedrals(traj,list_dih2)).flatten()[0]
  #main_dih,sub_dih2,sub_ang1=main_dih*180.0/math.pi,sub_dih2*180.0/math.pi,sub_ang1*180.0/math.pi
  main_dih,sub_ang1=main_dih*180.0/math.pi,sub_ang1*180.0/math.pi
  tg_class=tg_classtext(main_dih)
  #outfile.write("{} {:8.3f} {} {:8.3f} {:8.3f} {:8.3f} dim# NCCOdih {} {} {}\n".\
  outfile.write("{} {:8.3f} {} {:8.3f} {:8.3f} dim# NCCOdih {} {}\n".\
  #format(i,main_dih,tg_class,sub_dist1,sub_ang1,sub_dih2,text_sub_dist1,text_sub_ang1,text_sub_dih2))
  format(i,main_dih,tg_class,sub_dist1,sub_ang1,text_sub_dist1,text_sub_ang1))

outfile.close()

#  sub_dist1=np.minimum(sub_distc1,sub_distc2)
#  main_dih=(md.compute_dihedrals(traj,list_dih1)).flatten()[0]
#  sub_dih2=(md.compute_dihedrals(traj,list_dih2)).flatten()[0]
#  main_dih,sub_dih2,sub_ang1=main_dih*180.0/math.pi,sub_dih2*180.0/math.pi,sub_ang1*180.0/math.pi
  
#  tg_class=tg_classtext(main_dih)
#  outfile.write("{} {:8.3f} {} {:8.3f} {:8.3f} {:8.3f} dim# NCCOdih {} {} {}\n".\
#  format(i,main_dih,tg_class,sub_dist1,sub_ang1,sub_dih2,text_sub_dist1,text_sub_ang1,text_sub_dih2))
#outfile.close()
