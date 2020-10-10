#!/usr/bin/env python3
# nonbonding atomic pair potential energy function generator
# # transformed gromacs nonbonding potential table generator
# algorithm : this works for ChCl/urea switchable Cl force field only.
# put the cut xml file NB entries as input
# read atom types and constants
# generate pairwise nonbonding potential tables (atomtype combination : permutation)
# function of distance. Set interval : 0.1nm ~ 1.4nm?
# should be able to calculate for both OPLS, SAPT-FF.
# unit of energy : kJ/mol.
# better if the table supports separate columns of total, coulomb, LJ.

import math
import numpy as np
import argparse

#Desti = {'class' : -1, 'Aexch' : 0, 'Aelec' : 1, 'Aind' : 2, 'Adhf' : 3, 'Bexp' : 4,
#	 'C6' : 5, 'C8' : 6, 'C10' : 7, 'C12' : 8, 'Switch' : 9,
#	 'Aexadj' : 10, 'Aeladj' : 11, 'Ainadj' : 12, 'Adhadj' : 13}
ecp_const = 138.936 #electrostatic coupling constant: 138.936 kJ nm/(mol e^2)
parser = argparse.ArgumentParser(
  formatter_class=argparse.ArgumentDefaultsHelpFormatter, 
  description='generate nonbonding u(r) from openmm xml FF')
## args
parser.add_argument('-i', '--input', default='ff.xml', nargs='?', 
  help='input customNB partial xml file')

# parse xml input file
args = parser.parse_args()

fftype=input("Which FF is imported? OPLS? SAPT?(assuming non-pol version for charge)\n")

#if fftype == 'OPLS':
#  Desti = {'class' : -1, 'charge' : 14, 'sigma' : , 'epsilon' : }
#elif fftype == 'SAPT':
Desti = {'type' : -1, 'class' : -1, 'Aexch' : 0, 'Aelec' : 1, 'Aind' : 2, 'Adhf' : 3, 'Bexp' : 4,\
	 'C6' : 5, 'C8' : 6, 'C10' : 7, 'C12' : 8, 'Switch' : 9,\
	 'Aexadj' : 10, 'Aeladj' : 11, 'Ainadj' : 12, 'Adhadj' : 13,\
         'charge' : 14, 'sigma' : 15, 'epsilon' : 16}

f = open(args.input,'r')
lines = f.readlines()
data = [line.split() for line in lines]
natoms = len(data)
coefs=np.zeros((natoms,17))
aclass=[None]*natoms

for i in range(natoms):
  line=data[i]
  for word in line:
    if '=' in word:
      wsplit=word.split('=')
      prefix=wsplit[0]
      value=wsplit[1].split('"')[1]
      k=Desti[prefix]
      if k==-1:
        aclass[i]=value
      else:
        coefs[i][k]=np.float64(value)

fftable=np.zeros((651,4),dtype=np.float64)
r=np.array([x/1000.0 for x in range(100,1402,2)],dtype=np.float64)
fftable[:,0]=r

# Double loop of nb table.
for i in range(natoms):
  for j in range(i,natoms):
    #construct mutual parameters : OPLS or SAPT
    outname='nbur_'+aclass[i]+'_'+aclass[j]+'.dat'
    qi,qj=coefs[i][14],coefs[j][14]
    ucoul=ecp_const*qi*qj/r
    if fftype == 'OPLS':
      #OPLS uses geometric combination rule always
      si,sj=coefs[i][15],coefs[j][15]
      ei,ej=coefs[i][16],coefs[j][16]
      sij,eij=np.sqrt(si*sj),np.sqrt(ei*ej)
      g = -4.0*eij*((sij/r)**6)
      h = 4.0*eij*((sij/r)**12)
    elif fftype == 'SAPT':
      Aexch1,Aelec1,Aind1,Adhf1,Bexp1,C61,C81,C101,C121,Switch1,Aexadj1,Aeladj1,Ainadj1,Adhadj1= \
      coefs[i][0],coefs[i][1],coefs[i][2],coefs[i][3],coefs[i][4],coefs[i][5],coefs[i][6],\
      coefs[i][7],coefs[i][8],coefs[i][9],coefs[i][10],coefs[i][11],coefs[i][12],coefs[i][13]
      Aexch2,Aelec2,Aind2,Adhf2,Bexp2,C62,C82,C102,C122,Switch2,Aexadj2,Aeladj2,Ainadj2,Adhadj2= \
      coefs[j][0],coefs[j][1],coefs[j][2],coefs[j][3],coefs[j][4],coefs[j][5],coefs[j][6],\
      coefs[j][7],coefs[j][8],coefs[j][9],coefs[j][10],coefs[j][11],coefs[j][12],coefs[j][13]
      Aex=np.sqrt((Aexch1 + Aexadj1*(1-Switch1*Switch2))*(Aexch2 + Aexadj2*(1-Switch1*Switch2)))
      Ael=np.sqrt((Aelec1 + Aeladj1*(1-Switch1*Switch2))*(Aelec2 + Aeladj2*(1-Switch1*Switch2)))
      Ain=np.sqrt((Aind1 + Ainadj1*(1-Switch1*Switch2))*(Aind2 + Ainadj2*(1-Switch1*Switch2)))
      Adh=np.sqrt((Adhf1 + Adhadj1*(1-Switch1*Switch2))*(Adhf2 + Adhadj2*(1-Switch1*Switch2)))
      A=Aex-Ael-Ain-Adh
      B=(Bexp1+Bexp2)*Bexp1*Bexp2/(Bexp1**2 + Bexp2**2)
      C6,C8,C10,C12=np.sqrt(C61*C62),np.sqrt(C81*C82),np.sqrt(C101*C102),np.sqrt(C121*C122)

      Br=B*r
      exBr=np.exp(-Br,dtype=np.float64)
      f6 = 1 - exBr*(1 + Br * (1 + (1/2)*Br*(1 + (1/3)*Br*(1 + (1/4)*Br*(1 + (1/5)*Br*(1 + (1/6)*Br ) ) )  ) ) )
      #1-exBr(1+Br+1/2Br^2+1/6Br^3+1/24Br^4+ 1/120Br^5+1/720Br^6)
      f8 = f6 - exBr*( (1/5040)*(Br**7)*(1 + Br/8 ) )
      f10 = f8 - exBr*( (1/362880)*(Br**9)*(1 + Br/10 ) )
      f12 = f10 - exBr*( (1/39916800)*(Br**11)*(1 + Br/12) )
      #fp6=B*exBr*(1 + Br * (1 + (1/2)*Br*(1 + (1/3)*Br*(1 + (1/4)*Br*(1 + (1/5)*Br*(1 + (1/6)*Br ) ) )  ) ) ) \
      #-exBr*(B+(B**2)*r+(1/2)*(B**3)*(r**2)+(1/6)*(B**4)*(r**3)+(1/24)*(B**5)*(r**4)+(1/120)*(B**6)*(r**5))
      #fp8=fp6 + B*exBr*((1/5040)*(Br**7)*(1 + Br/8 ))-exBr*( (1/720)*(B**7)*(r**6)+(1/5040)*(B**8)*(r**7) )
      #fp10=fp8 + B*exBr*( (1/362880)*(Br**9)*(1 + Br/10 ) ) - exBr*( (1/40320)*(B**9)*(r**8)+(1/362880)*(B**10)*(r**9) )
      #fp12=fp10 + B*exBr*( (1/39916800)*(Br**11)*(1 + Br/12) ) - exBr*( (1/3628800)*(B**11)*(r**10)+(1/39916800)*(B**12)*(r**11))
      g = -f6*C6/(r**6)-f8*C8/(r**8)-f10*C10/(r**10)-f12*C12/(r**12)
      #gprime = (-6)*f6*C6/(r**7)-8*f8*C8/(r**9)-10*f10*C10/(r**11)-12*f12*C12/(r**13) + fp6*C6/(r**6)+fp8*C8/(r**8)+fp10*C10/(r**10)+fp12*C12/(r**12)
      #h,hprime=A*exBr,A*B*exBr
      h = A*exBr
      #fftable[:,3],fftable[:,4],fftable[:,5],fftable[:,6]=g,gprime,h,hprime
    unocoul=g+h
    utot=ucoul+unocoul
    fftable[:,1],fftable[:,2],fftable[:,3]=utot,ucoul,unocoul
    #print table file
    fffile=open(outname,'w')
    for line in fftable:
      #prevent scale-dominating data points: only accept within +- 1000 kJ/mol(totE)
      if line[1]<=1000 and line[1]>=-1000:
        fffile.write('{:12.10e} {:12.10e} {:12.10e} {:12.10e}\n'.\
        format(line[0],line[1],line[2],line[3]))
    fffile.close()

#compute tablated potential
# <CustomNonbondedForce energy="A*exBr - f6*C6/(r^6) - f8*C8/(r^8) - f10*C10/(r^10) - f12*C12/(r^12);
#    A=Aex-Ael-Ain-Adh;
#    Aex=sqrt((Aexch1 + Aexadj1*(1-Switch1*Switch2))*(Aexch2 + Aexadj2*(1-Switch1*Switch2)));
#    Ael=sqrt((Aelec1 + Aeladj1*(1-Switch1*Switch2))*(Aelec2 + Aeladj2*(1-Switch1*Switch2)));
#    Ain=sqrt((Aind1 + Ainadj1*(1-Switch1*Switch2))*(Aind2 + Ainadj2*(1-Switch1*Switch2)));
#    Adh=sqrt((Adhf1 + Adhadj1*(1-Switch1*Switch2))*(Adhf2 + Adhadj2*(1-Switch1*Switch2)));
#    f12 = f10 - exBr*( (1/39916800)*(Br^11)*(1 + Br/12) );
#    f10 = f8 - exBr*( (1/362880)*(Br^9)*(1 + Br/10 ) );
#    f8 = f6 - exBr*( (1/5040)*(Br^7)*(1 + Br/8 ) );
#    f6 = 1 - exBr*(1 + Br * (1 + (1/2)*Br*(1 + (1/3)*Br*(1 + (1/4)*Br*(1 + (1/5)*Br*(1 + (1/6)*Br ) ) )  ) ) );
#    exBr = exp(-Br);
#    Br = B*r;
#    B=(Bexp1+Bexp2)*Bexp1*Bexp2/(Bexp1^2 + Bexp2^2);
#    C6=sqrt(C61*C62); C8=sqrt(C81*C82); C10=sqrt(C101*C102); C12=sqrt(C121*C122)"
### GROMACS tabulated potential rule ##
# r : grid of 0.002nm (for single-precision). 1501 grid points (toward 3.000 nm)
# 7 columns: consist of r,f,f',g,g',h,h'  . Also, values are defined to be 0 when r<=0.04nm. in C language, number format is %12.10e
# f,g,r functions are defined from : V(r) = (qiqj/4pie0) f(r) + Cg(r) + Ah(r)  (elst, dispersion, short range repulsion)
# parameters A,C can be specified in topology file. if using custom potentials, these parameters should be set to 1. 
# for LJ 6-12 potential,  f = 1/r, f' = 1/r^2, g=-1/r^6, g'=-6/r^7, h=1/r^12, 12/r^13
#for sapt FF: g= - f6*C6/(r^6) - f8*C8/(r^8) - f10*C10/(r^10) - f12*C12/(r^12), g' = -6*f6*C6/(r^7) - 8*f8*C8/(r^9) - 10*f10*C10/(r^11) - 12*f12*C12/(r^13)
#for sapt FF : h = A*exBr, h' = A*B*exBr

#note : Coulomb potential energy
# u_coulomb = qiqj e^2 /(4pie0)rij
# 1/(4pie0) : Coulomb const., 8.98755 x 10^9 N m2 C-2
# e^2 elementary charge^2 : 2.56697 x 10^-38 C2
# qiqj can be in elementary charge unit. 
# For rij in nm (10^-9 m), the coefficient is : 2.30708 x 10 * 10^-38 * 10^9 / 10^-9 = 
# 2.30708 x 10^-19 J .
# in kJ/mol : 2.30708 x 10^-19 * 6.02214 x 10^23 / 1000 = 138.936 kJ/mol
# It is consistent with published 'electrostatic coupling constant' 138.936 kJ nm/(mol e^2).   
