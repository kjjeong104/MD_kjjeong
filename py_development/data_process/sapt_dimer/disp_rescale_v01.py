#!/home/kjeong23/softwares/bin/python3.4
#SAPT FF dispersion coefficient rescaling program
#algorithm : read SAPT data : coord and energy.
#recalculate (C6-C12)FF energy from coordinate.
#try to fit with C12=0. For the formula, see SAPT-FF papers.(be careful about damping factors)
#should read ab-initio dispersion energy as well.
#ab-initio dispersion energy: E2disp + E2dispexch
#FF dispersion energy: -(sigma_6,8,10,12)(sigma_ij)f(Bij,rij)Cnij/rnij
#Tang-Toennies damping function fn(B,r) = 1-exp(-Br)(sigma_m=0~n)((Br)^m/m!)
#first test if the FF energy calculation is done correctly.

import math
import sys
import numpy
from scipy.optimize import curve_fit
from scipy.special import factorial

def damp(b1,b2,r,n):
  l=(b1+b2)*(b1*b2)/(b1**2+b2**2)
  s=1.0
  for i in range(1,n+1):
    s+=((l*r)**i)/float(factorial(i))
  damp=1.0-numpy.exp(-l*r)*s
  return damp

def ff_ecomp_disp_calc(crd1,crd2,B1,B2,Cn1,Cn2): #comprehensive calc of FF disp energy
  # for all snapshots. Gives array of FF disp energies
  #crd i : all snapshots of monomer i crd. 3-dim array [snapshot][atom][xyz]
  nstep=crd1.shape[0] # number of snapshots
  natom1,natom2=crd1.shape[1],crd2.shape[1]
  print(natom1,natom2)
  Cn_cross=numpy.empty((Cn1.shape[0],Cn2.shape[0],4))
  for i in range(Cn1.shape[0]):
    for j in range(Cn2.shape[0]):
      Cn_cross[i][j]=numpy.sqrt(Cn1[i]*Cn2[j])

  ff_ecomp_disp=numpy.zeros(nstep)
  for n in range(nstep):
    s=0.0 #total intermolecular ff dispersion energy of 1 snapshot
    for i in range(natom1):
      for j in range(natom2):
        rij=crd1[n][i]-crd2[n][j]
        dist2=numpy.dot(rij,rij)
        dist=numpy.sqrt(dist2)
        for k in range(1,5): #C6,C8,C10,C12 contrib
          order=6+(k-1)*2
          term=damp(B1[i],B2[j],dist,order)*Cn_cross[i][j][k-1]/dist2**(order/2)
          s=s-term
    ff_ecomp_disp[n]=s
  return ff_ecomp_disp

def main():
  fffile=open(sys.argv[1],'r') #input parameters
  bohrfile=open(sys.argv[2],'r') #sapt calc results
  outfile= open(sys.argv[3],'w') #sapt vs ff comparison plot

  aname1,aname2=[],[] #monomer1 and monomer2
  B1,B2=numpy.empty(0),numpy.empty(0)
  Cn1,Cn2=numpy.empty((0,4)),numpy.empty((0,4))

  #read parameters : exp, disp.
  mflag,pflag=0,0 #monomer flag, parameter type flag
  while True:
    ffline=fffile.readline()
    if ffline=='':#EOF
      break
    elif 'mon1' in ffline:#monomer 1
      mflag=1
    elif 'mon2' in ffline:#monomer 2
      mflag=2
    elif 'exponents' in ffline:
      pflag=1
    elif 'dispersion' in ffline:
      pflag=2
    else:
      ffsplit=ffline.split()
      if pflag==1:
        if mflag==1:
          aname1.append(ffsplit[0])
          B1=numpy.append(B1,float(ffsplit[1]))
        elif mflag==2:
          aname2.append(ffsplit[0])
          B2=numpy.append(B2,float(ffsplit[1]))
      elif pflag==2:
        newrow=numpy.array([float(ffsplit[1]),float(ffsplit[2]),float(ffsplit[3]),float(ffsplit[4])])
        if mflag==1:
          Cn1=numpy.vstack((Cn1,newrow))
        elif mflag==2:
          Cn2=numpy.vstack((Cn2,newrow))
  ##print(B1,B2)
  ##print(Cn1,Cn2)

  #read bohrfile : get coordinate, energy value
  sindex,natom1,natom2=0,0,0
  dispe=0.0
  crd1,crd2=numpy.empty((0,0,3)),numpy.empty((0,0,3))
  sapt_ecomp_disp=numpy.empty(0)
  while True:
    bohrline=bohrfile.readline()
    if bohrline=='':#EOF
      break
    bsplit=bohrline.split()
    if len(bsplit)==1: #num of atoms
      if natom1==0 and natom2==0: #both undefined
        natom1=int(bsplit[0])
        crdtemp=numpy.empty((natom1,3))
        for i in range(natom1): #read 1 monomer traj (mon1)
          bohrline=bohrfile.readline()
          bsplit=bohrline.split()
          newrow=numpy.array([float(bsplit[1]),float(bsplit[2]),float(bsplit[3])])
          crdtemp[i]=newrow
        crdtemp=numpy.array([crdtemp]) #make 3-d array
        if sindex==0:
          crd1=crdtemp
        else:
          crd1=numpy.append(crd1,crdtemp,axis=0)

      elif natom1!=0 and natom2==0:
        natom2=int(bsplit[0])
        crdtemp=numpy.empty((natom2,3))
        for j in range(natom2): #read 1 monomer traj (mon2)
          bohrline=bohrfile.readline()
          bsplit=bohrline.split()
          newrow=numpy.array([float(bsplit[1]),float(bsplit[2]),float(bsplit[3])])
          crdtemp[j]=newrow
        crdtemp=numpy.array([crdtemp]) #make 3-d array
        if sindex==0:
          crd2=crdtemp
        else:
          crd2=numpy.append(crd2,crdtemp,axis=0)

    elif len(bsplit)>1:
      if bsplit[0]=='E2disp':
        dispe+=float(bsplit[1])/1000.0
      elif bsplit[0]=='E2disp-exch':
        dispe+=float(bsplit[1])/1000.0
        sapt_ecomp_disp=numpy.append(sapt_ecomp_disp,dispe)
        dispe=0.0
      elif 'dhf' in bohrline: #call for end of a step
        sindex+=1
        natom1,natom2=0,0

  ##print(sapt_ecomp_disp)
  #calculate
  ff_ecomp_disp=ff_ecomp_disp_calc(crd1,crd2,B1,B2,Cn1,Cn2)

  #writing section (sapt vs ff)
  for i in range(len(ff_ecomp_disp)):
    outfile.write('{:15.8f} {:15.8f}\n'.format(sapt_ecomp_disp[i],ff_ecomp_disp[i]))

  bohrfile.close()
  fffile.close()

if __name__ == "__main__": main()
