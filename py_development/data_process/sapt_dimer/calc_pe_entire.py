from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import numpy as np

# load fake pdb file. Later, new modeller will be defined, and positions will be updated
# do single point energy calculation for mm-relaxed geometries, to calculate residual PE

temperature=300*kelvin
pdbtemp="qm_HCCO_0.pdb"
ffin=sys.argv[1] #newly fit NCCO parameters. zeroed HCCO. or use arguments
pdb = PDBFile(pdbtemp)
strdir = ''

integ_md = DrudeLangevinIntegrator(temperature, 1/picosecond, 1*kelvin, 1/picosecond, 0.001*picoseconds)
integ_md.setMaxDrudeDistance(0.02)  # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)
pdb.topology.loadBondDefinitions('sapt_residues_choline.xml')
pdb.topology.createStandardBonds();

modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField(ffin)
modeller.addExtraParticles(forcefield)
#PDBFile.writeFile(modeller.topology, modeller.positions, open('init.pdb', 'w'))

system = forcefield.createSystem(modeller.topology, constraints=None, rigidWater=True)
nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]
customNonbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomNonbondedForce][0]
drudeForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == DrudeForce][0]
nbondedForce.setNonbondedMethod(NonbondedForce.NoCutoff)
customNonbondedForce.setNonbondedMethod(min(nbondedForce.getNonbondedMethod(),NonbondedForce.NoCutoff))

for i in range(system.getNumForces()):
    f = system.getForce(i)
    type(f)
    f.setForceGroup(i)

totmass = 0.*dalton
for i in range(system.getNumParticles()):
    totmass += system.getParticleMass(i)
simmd = Simulation(modeller.topology, system, integ_md)

for a in range(0,360,3):
  uppdb = PDBFile("mm_NCCO_"+str(a)+".pdb")
  modeller2=Modeller(modeller.topology,uppdb.positions)
  simmd.context.setPositions(modeller2.positions)
  state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
  #position = state.getPositions()
  print('NCCO dih :',a,' PE(total) : '+str(state.getPotentialEnergy()))
  for i in range(system.getNumForces()):
    f = system.getForce(i)
    print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()))


#for i in range(system.getNumForces()):
#    f = system.getForce(i)
#    print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()))

#PDBFile.writeFile(simmd.topology, position, open(strdir+'beforemin.pdb', 'w'))
#print('Wrote initial positions')

print('Done!')

exit()
