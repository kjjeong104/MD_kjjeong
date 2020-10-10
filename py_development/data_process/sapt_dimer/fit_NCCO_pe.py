from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import numpy as np
from scipy.optimize import curve_fit

# load fake pdb file. Later, new modeller will be defined, and positions will be updated
# do single point energy calculation for mm-relaxed geometries, to calculate residual PE

temperature=300*kelvin
pdbtemp="qm_HCCO_0.pdb"
ffin="ff_onlyNCCO.xml" #ff with only HCCO dihedral
pdb = PDBFile(pdbtemp)
strdir = ''

infile=sys.argv[1]
odata=numpy.loadtxt(infile)
x=odata[:,0]
targety=odata[:,1]

integ_md = LangevinIntegrator(temperature, 1/picosecond, 0.001*picoseconds)
pdb.topology.loadBondDefinitions('sapt_residues_choline.xml')
pdb.topology.createStandardBonds();

modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField(ffin)
modeller.addExtraParticles(forcefield)

system = forcefield.createSystem(modeller.topology, constraints=None, rigidWater=True)
torsion = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == RBTorsionForce][0]

for i in range(system.getNumForces()):
    f = system.getForce(i)
    type(f)
    f.setForceGroup(i)

totmass = 0.*dalton
for i in range(system.getNumParticles()):
    totmass += system.getParticleMass(i)
simmd = Simulation(modeller.topology, system, integ_md)

#for a in range(0,360,3):
#  uppdb = PDBFile("mm_NCCO_"+str(a)+".pdb")
#  modeller2=Modeller(modeller.topology,uppdb.positions)
#  simmd.context.setPositions(modeller2.positions)
#  state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
#  #position = state.getPositions()
#  print('HCCO dih :',a,' PE(residual) : '+str(state.getPotentialEnergy()))

p1,p2,p3,p4=0,13,16,19
def function(x,c0,c1,c2,c3,c4,c5):
  torsion.setTorsionParameters(0,p1,p2,p3,p4,c0,c1,c2,c3,c4,c5)
  torsion.updateParametersInContext(simmd.context)
  y=numpy.zeros(len(x))
  i=0
  for entry in x:
    uppdb= PDBFile("mm_NCCO_"+str(int(entry))+".pdb")
    modeller2=Modeller(modeller.topology,uppdb.positions)
    simmd.context.setPositions(modeller2.positions)
    state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
    y[i]=state.getPotentialEnergy()/kilojoules*mole #PES as function of x, using parameters c0~c5
    i+=1

  return y

new_coeffs,cov=curve_fit(function,x,targety,method='lm')
print(new_coeffs)

y=function(x,new_coeffs[0],new_coeffs[1],new_coeffs[2],new_coeffs[3],new_coeffs[4],new_coeffs[5])
#print(y)
for i in range(len(x)):
  print(x[i],y[i],targety[i])

print('Done!')

exit()
