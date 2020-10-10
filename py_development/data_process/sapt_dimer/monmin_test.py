from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import numpy as np

def angle(p):
    a,b,c=p[0],p[1],p[2]
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

def dihedral(p):
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2( y, x ))

temperature=300*kelvin
pdbin=sys.argv[1]
ffin=sys.argv[2]
pdbout=sys.argv[3]
pdb = PDBFile(pdbin)
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
print('nbMethod : ', customNonbondedForce.getNonbondedMethod())

for i in range(system.getNumForces()):
    f = system.getForce(i)
    type(f)
    f.setForceGroup(i)

totmass = 0.*dalton
for i in range(system.getNumParticles()):
    totmass += system.getParticleMass(i)
simmd = Simulation(modeller.topology, system, integ_md)
simmd.context.setPositions(modeller.positions)
simmd.context.reinitialize()
simmd.context.setPositions(modeller.positions)

print('Minimizing...')
state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
position = state.getPositions()
nppos=state.getPositions(asNumpy=True)
print(str(state.getKineticEnergy()))
print(str(state.getPotentialEnergy()))

a1atoms=nppos[[14,13,15],:]/(1.0*nanometer)
a2atoms=nppos[[18,17,19],:]/(1.0*nanometer)
a3atoms=nppos[[0,13,16],:]/(1.0*nanometer)
d1atoms=nppos[[14,13,16,19],:]/(1.0*nanometer)
d2atoms=nppos[[0,13,16,19],:]/(1.0*nanometer)
print('chain Bond angle HCH (beforemin) : ',angle(a1atoms),' ',angle(a2atoms))
print('Bond angle NCC (beforemin) : ',angle(a3atoms))
print('Dihedral angle HCCO (beforemin) : ',dihedral(d1atoms))
print('Dihedral angle NCCO (beforemin) : ',dihedral(d2atoms))

for i in range(system.getNumForces()):
    f = system.getForce(i)
    print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()))

#PDBFile.writeFile(simmd.topology, position, open(strdir+'beforemin.pdb', 'w'))
print('Wrote initial positions')
simmd.minimizeEnergy(maxIterations=2000)

print('Minimization finished !')
state = simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
print(numpy.max(state.getVelocities()*picosecond/nanometer))
print(str(state.getKineticEnergy()))
print(str(state.getPotentialEnergy()))
for i in range(system.getNumForces()):
    f = system.getForce(i)
    print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()))

position = state.getPositions()
#print(position)
nppos=state.getPositions(asNumpy=True)
a1atoms=nppos[[14,13,15],:]/(1.0*nanometer)
a2atoms=nppos[[18,17,19],:]/(1.0*nanometer)
a3atoms=nppos[[0,13,16],:]/(1.0*nanometer)
d1atoms=nppos[[14,13,16,19],:]/(1.0*nanometer)
d2atoms=nppos[[0,13,16,19],:]/(1.0*nanometer)
print('chain Bond angle HCH (aftermin) : ',angle(a1atoms),' ',angle(a2atoms))
print('Bond angle NCC (beforemin) : ',angle(a3atoms))
print('Dihedral angle HCCO (aftermin) : ',dihedral(d1atoms))
print('Dihedral angle NCCO (aftermin) : ',dihedral(d2atoms))

PDBFile.writeFile(simmd.topology, position, open(strdir+pdbout, 'w'))

#*************************************
#*************************************************
# ChangYun created the DrudeDataReporter class, we need to pull this into this OpenMM install if we want to use it
#*************************************************
#simmd.reporters.append(DrudeDataReporter(strdir+'md_nvt.log', 1000, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, langevin=True, density=False,speed=True))
#simmd.reporters.append(DrudeDataReporter(strdir+'md_nvt_temp.log', 10000, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, langevin=True, drudeTemperature=True,density=False,speed=True))
#simmd.reporters[2].report(simmd,state)

#for i in range(simmd.system.getNumForces()):
#    if type(simmd.system.getForce(i)) == MonteCarloBarostat:
#        simmd.system.removeForce(i)
print('Done!')

exit()
