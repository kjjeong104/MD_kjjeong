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
#v01 : only showing current disp energy
#v02 : fitting with zeroing C12 term of an atom
#v021 : uniform scale-up for C6,C8,C10
#v022 : merge v02 and v022 functions. Also, as v01 does, show updated plot for (sapt vs ff fit)

import math
import sys
import numpy
import timeit
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
  start_time=timeit.default_timer()

  aname1,aname2=[],[] #monomer1 and monomer2
  B1,B2=numpy.empty(0),numpy.empty(0)
  Cn1,Cn2=numpy.empty((0,4)),numpy.empty((0,4))

  #read parameters : exp, disp.
  mflag,pflag=0,0 #monomer flag, parameter type flag
  pre_index1,pre_index2=-1,-1
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
          pre_index1+=1
        elif mflag==2:
          Cn2=numpy.vstack((Cn2,newrow))
          pre_index2+=1
        if len(ffsplit)==6:
          if ffsplit[5]=='C12_zero':
            newdisp_mon=mflag
            if mflag==1:
              newdisp_atom=pre_index1
            elif mflag==2:
              newdisp_atom=pre_index2
            guessc6810=numpy.array([float(ffsplit[1]),float(ffsplit[2]),float(ffsplit[3])])
  ##print(B1,B2)
  ##print(Cn1,Cn2)
  print(newdisp_mon,newdisp_atom)

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
          crd1=numpy.copy(crdtemp)
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
          crd2=numpy.copy(crdtemp)
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
  ff_ecomp_disp_ori=ff_ecomp_disp_calc(crd1,crd2,B1,B2,Cn1,Cn2)

  x=numpy.arange(len(ff_ecomp_disp_ori))
  #function to opt new coef
  #target for new fitting : monomer (newdisp_mon), atom index (start from 0: newdisp_atom)
  def newdisp_rescale(x,scale): #new disp calculation for 1 snapshot. C12 discarded.
    #crdtemp1,crdtemp2=crd1[x],crd2[x]
    nstep=len(x)
    natom1,natom2=crd1.shape[1],crd2.shape[1]
    Cn_cross1=numpy.empty((Cn1.shape[0],Cn2.shape[0],4))
    Cn1temp,Cn2temp=numpy.copy(Cn1),numpy.copy(Cn2)
    #overwrite new Cn parameters for fitting target.
    if newdisp_mon==1:
      Cn1temp[newdisp_atom]=Cn1[newdisp_atom]*scale
      Cn1temp[newdisp_atom][3]=0.0
    elif newdisp_mon==2:
      Cn2temp[newdisp_atom]=Cn2[newdisp_atom]*scale
      Cn2temp[newdisp_atom][3]=0.0
    #construct combined Cn
    for i in range(Cn1temp.shape[0]):
      for j in range(Cn2temp.shape[0]):
        Cn_cross1[i][j]=numpy.sqrt(Cn1temp[i]*Cn2temp[j])

    y=numpy.zeros(nstep)
    n=0
    for entry in x:
      s=0.0 #total intermolecular ff dispersion energy of 1 snapshot
      for i in range(natom1):
        for j in range(natom2):
          rij=crd1[entry][i]-crd2[entry][j]
          dist2=numpy.dot(rij,rij)
          dist=numpy.sqrt(dist2)
          term6=damp(B1[i],B2[j],dist,6)*Cn_cross1[i][j][0]/dist2**(6/2)
          term8=damp(B1[i],B2[j],dist,8)*Cn_cross1[i][j][1]/dist2**(8/2)
          term10=damp(B1[i],B2[j],dist,10)*Cn_cross1[i][j][2]/dist2**(10/2)
          s=s-term6-term8-term10
      y[n]=s
      n+=1

    elapsed=timeit.default_timer() - start_time
    print('iterating.. time {:15.8f}'.format(elapsed))
    print('candidate scale: {:15.8f}'.format(scale))
    return y
  def newdisp_refit(x,c6,c8,c10): #new disp calculation for 1 snapshot. C12 discarded.
    #crdtemp1,crdtemp2=crd1[x],crd2[x]
    nstep=len(x)
    natom1,natom2=crd1.shape[1],crd2.shape[1]
    Cn_cross1=numpy.empty((Cn1.shape[0],Cn2.shape[0],4))
    Cn1temp,Cn2temp=numpy.copy(Cn1),numpy.copy(Cn2)
    #overwrite new Cn parameters for fitting target.
    if newdisp_mon==1:
      Cn1temp[newdisp_atom]=numpy.array([c6,c8,c10,0.0])
    elif newdisp_mon==2:
      Cn2temp[newdisp_atom]=numpy.array([c6,c8,c10,0.0])
    #construct combined Cn
    for i in range(Cn1temp.shape[0]):
      for j in range(Cn2temp.shape[0]):
        Cn_cross1[i][j]=numpy.sqrt(Cn1temp[i]*Cn2temp[j])

    y=numpy.zeros(nstep)
    n=0
    for entry in x:
      s=0.0 #total intermolecular ff dispersion energy of 1 snapshot
      for i in range(natom1):
        for j in range(natom2):
          rij=crd1[entry][i]-crd2[entry][j]
          dist2=numpy.dot(rij,rij)
          dist=numpy.sqrt(dist2)
          term6=damp(B1[i],B2[j],dist,6)*Cn_cross1[i][j][0]/dist2**(6/2)
          term8=damp(B1[i],B2[j],dist,8)*Cn_cross1[i][j][1]/dist2**(8/2)
          term10=damp(B1[i],B2[j],dist,10)*Cn_cross1[i][j][2]/dist2**(10/2)
          s=s-term6-term8-term10
      y[n]=s
      n+=1

    elapsed=timeit.default_timer() - start_time
    print('iterating.. time {:15.8f}'.format(elapsed))
    print('candidate c6,c8,c10 coefs: {:15.8f} {:15.8f} {:15.8f}'.format(c6,c8,c10))
    return y

  #new_scale,cov=curve_fit(newdisp,x,ff_ecomp_disp,p0=1.0,method='lm')
  new_scale,cov=curve_fit(newdisp_rescale,x,sapt_ecomp_disp,p0=1.0,method='lm')
  new_coeffs,cov=curve_fit(newdisp_refit,x,sapt_ecomp_disp,p0=guessc6810,method='lm')
  #modify coefs, and calculate ecomp_disp again
  Cn1rescale,Cn2rescale=numpy.copy(Cn1),numpy.copy(Cn2)
  Cn1refit,Cn2refit=numpy.copy(Cn1),numpy.copy(Cn2)
  if newdisp_mon==1:
    Cn1rescale[newdisp_atom]=Cn1rescale[newdisp_atom]*new_scale
    Cn1rescale[newdisp_atom][3]=0.0
    Cn1refit[newdisp_atom]=numpy.array([new_coeffs[0],new_coeffs[1],new_coeffs[2],0.0])
  elif newdisp_mon==2:
    Cn2rescale[newdisp_atom]=Cn2rescale[newdisp_atom]*new_scale
    Cn2rescale[newdisp_atom][3]=0.0  
    Cn2refit[newdisp_atom]=numpy.array([new_coeffs[0],new_coeffs[1],new_coeffs[2],0.0])
  ff_ecomp_disp_rescale=ff_ecomp_disp_calc(crd1,crd2,B1,B2,Cn1rescale,Cn2rescale) 
  ff_ecomp_disp_refit=ff_ecomp_disp_calc(crd1,crd2,B1,B2,Cn1refit,Cn2refit) 

  #writing section (sapt vs ff)
  for i in range(len(ff_ecomp_disp_ori)):
    outfile.write('{:15.8f} {:15.8f} {:15.8f} {:15.8f}\n'.format(sapt_ecomp_disp[i],ff_ecomp_disp_ori[i],ff_ecomp_disp_rescale[i],ff_ecomp_disp_refit[i]))
  print('new scale: ',new_scale)
  print('scaled coefs: ',guessc6810*new_scale)
  print('if refitted, new coefs: ',new_coeffs)

  bohrfile.close()
  fffile.close()
  outfile.close()

if __name__ == "__main__": main()
