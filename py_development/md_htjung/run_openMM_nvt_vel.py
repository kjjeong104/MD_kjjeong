from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime

temperature=323.15*kelvin

pdb = PDBFile('nvt_init.pdb')
strdir = ''
cpu_gpu= 'CPU'

integ_md = DrudeLangevinIntegrator(temperature, 1/picosecond, 1*kelvin, 1/picosecond, 0.001*picoseconds)
integ_md.setMaxDrudeDistance(0.02)  # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)

pdb.topology.loadBondDefinitions('sapt_residues.xml')
pdb.topology.createStandardBonds();

modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('sapt.xml')
modeller.addExtraParticles(forcefield)
#PDBFile.writeFile(modeller.topology, modeller.positions, open('init.pdb', 'w'))

system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=None, rigidWater=True)
nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]
customNonbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomNonbondedForce][0]
custombond = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomBondForce][0]
drudeForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == DrudeForce][0]
nbondedForce.setNonbondedMethod(NonbondedForce.PME)
customNonbondedForce.setNonbondedMethod(min(nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic))
print('nbMethod : ', customNonbondedForce.getNonbondedMethod())

for i in range(system.getNumForces()):
    f = system.getForce(i)
    type(f)
    f.setForceGroup(i)

totmass = 0.*dalton
for i in range(system.getNumParticles()):
    totmass += system.getParticleMass(i)

#system.addForce(AndersenThermostat(300*kelvin, 10/picosecond))

if 'CPU' in cpu_gpu:
  platform = Platform.getPlatformByName('CPU')
  simmd = Simulation(modeller.topology, system, integ_md, platform)
elif 'CL' in cpu_gpu:  
  platform = Platform.getPlatformByName('OpenCL')
  properties = {'OpenCLPrecision': 'mixed'}
  simmd = Simulation(modeller.topology, system, integ_md, platform, properties)
else:
  platform = Platform.getPlatformByName('CUDA')
  properties = {'CudaPrecision': 'mixed'}
  simmd = Simulation(modeller.topology, system, integ_md, platform, properties)

simmd.context.setPositions(modeller.positions)
#simmd.context.setPositions(pdbpos.positions)

platform =simmd.context.getPlatform()
platformname = platform.getName();
print(platformname)


#*************************************************************************************
#   This is section to create exclusions for TFSI nonbonded interactions, and update
#   Screened Drude interactions.  1-5 non-Coulomb interaction are accounted for
#   using CustomBondForce
#*************************************************************************************
print('Creating Exclusions for TFSI')

# map from global particle index to drudeforce object index
particleMap = {}
for i in range(drudeForce.getNumParticles()):
  particleMap[drudeForce.getParticleParameters(i)[0]] = i

# can't add duplicate ScreenedPairs, so store what we already have
flagexceptions = {}
for i in range(nbondedForce.getNumExceptions()):
  (particle1, particle2, charge, sigma, epsilon) = nbondedForce.getExceptionParameters(i)
  string1=str(particle1)+"_"+str(particle2)
  string2=str(particle2)+"_"+str(particle1)
  flagexceptions[string1]=1
  flagexceptions[string2]=1

# can't add duplicate customNonbonded exclusions, so store what we already have
flagexclusions = {}
for i in range(customNonbondedForce.getNumExclusions()):
  (particle1, particle2) = customNonbondedForce.getExclusionParticles(i)
  string1=str(particle1)+"_"+str(particle2)
  string2=str(particle2)+"_"+str(particle1)
  flagexclusions[string1]=1
  flagexclusions[string2]=1

# add exclusions for all atom pairs on TFSI residues, and when a drude pair is 
# excluded add a corresponding screened thole interaction in its place
for res in simmd.topology.residues():
  if res.name == 'Tf2N':
    for i in range(len(res._atoms)-1):
        for j in range(i+1,len(res._atoms)):
            (indi,indj) = (res._atoms[i].index, res._atoms[j].index)
	    # here it doesn't matter if we already have this, since we pass the "True" flag
            nbondedForce.addException(indi,indj,0,1,0,True)
	    # make sure we don't already exlude this customnonbond
            string1=str(indi)+"_"+str(indj)
            string2=str(indj)+"_"+str(indi)
            if string1 in flagexclusions and string2 in flagexclusions:
               continue
            else:
               customNonbondedForce.addExclusion(indi,indj)
            print(res.index,i,j)
            # add thole if we're excluding two drudes
            if indi in particleMap and indj in particleMap:
               # make sure we don't already have this screened pair
               if string1 in flagexceptions or string2 in flagexceptions:
                 continue
               else:
                 drudei = particleMap[indi]
                 drudej = particleMap[indj]
                 drudeForce.addScreenedPair(drudei, drudej, 2.0)


#************************************* end of TFSI section **************************
simmd.context.reinitialize()
simmd.context.setPositions(modeller.positions)

# load restart files
#simmd.loadCheckpoint('md_nvt.chk')

print('Starting Production NVT Simulation...')
t1 = datetime.now()
state = simmd.context.getState(getEnergy=True,getForces=False,getVelocities=False,getPositions=True)

# write initial pdb with drude oscillators
position = state.getPositions()
simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
PDBFile.writeFile(simmd.topology, position, open(strdir+'start_drudes.pdb', 'w'))

# initial energies
print(str(state.getKineticEnergy()))
print(str(state.getPotentialEnergy()))
for j in range(system.getNumForces()):
   f = system.getForce(j)
   print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))


simmd.reporters = []
simmd.reporters.append(DCDReporter(strdir+'md_nvt_pos.dcd',5,enforcePeriodicBox=True,writeVel=False))
simmd.reporters.append(DCDReporter(strdir+'md_nvt_vel.dcd',5,enforcePeriodicBox=False,writeVel=True))

#*************************************************
# ChangYun created the DrudeDataReporter class, we need to pull this into this OpenMM install if we want to use it
#*************************************************
#simmd.reporters.append(DrudeDataReporter(strdir+'md_nvt.log', 1000, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, langevin=True, density=False,speed=True))
#simmd.reporters.append(DrudeDataReporter(strdir+'md_nvt_temp.log', 10000, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, langevin=True, drudeTemperature=True,density=False,speed=True))
simmd.reporters.append(CheckpointReporter(strdir+'md_nvt.chk', 10000))
simmd.reporters[1].report(simmd,state)
#simmd.reporters[2].report(simmd,state)

#for i in range(simmd.system.getNumForces()):
#    if type(simmd.system.getForce(i)) == MonteCarloBarostat:
#        simmd.system.removeForce(i)

print('Simulating...')

for i in range(1,11):
    simmd.step(5)
    #state = simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,enforcePeriodicBox=True)
    #position = state.getPositions(asNumpy=True)
    #pos = position._value # nm
    #velocity = state.getVelocities(asNumpy=True)
    #vel = velocity.in_units_of(meter/second)._value # m/s
    #print(pos[0])
    #print(vel[0])
    #print(i,strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    #print(i,datetime.now())
    #state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True)
    #print(str(state.getKineticEnergy()))
    #print(str(state.getPotentialEnergy()))
    #for j in range(system.getNumForces()):
    #    f = system.getForce(j)
    #    print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))

t2 = datetime.now()
t3 = t2 - t1
print('simulation took', t3.seconds,'seconds')
print('Simulating...')
state = simmd.context.getState(getEnergy=True,getForces=True,getPositions=True,enforcePeriodicBox=True)
position = state.getPositions()
simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
PDBFile.writeFile(simmd.topology, position, open(strdir+'md_npt.pdb', 'w'))

print('Done!')

exit()
