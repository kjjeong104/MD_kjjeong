# cohesive energy calculation
# should use NVT simulation dcd file 
# and keep in mind that you need to modify pdb file (see below)
# 

from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import sys

############## variable settings #######################################################
temperature=298.15*kelvin
#pressure = Vec3(1.0,1.0,1.0)*atmosphere
#pressure = 1.0*atmosphere
#pressure = 0.986923*atmosphere
#barofreq = 100
# should run followings in command before using nvt_innit.pdb:
# sed -i '/D  /d' nvt_init.pdb
# sed -i '/CONECT/d' nvt_init.pdb
# sed -i -e 's/Cho A/CholA/g' nvt_init.pdb
# sed -i -e 's/Tf2 B/Tf2NB/g' nvt_init.pdb
init_pdb=str(sys.argv[1])  # first argument for pdb filename without drude particle
strdir = '' # path
def_bonds_xml='sapt_residues.xml' # definition xml file for bonds of residues
def_forces_xml='sapt_cohesive.xml' # definition xml file for bonding and non-bonding forces 
cpu_gpu_opencl='cpu'

cohesive_dcd=str(sys.argv[2]) # second argument for dcd trajectory file
cohesive_top=str(sys.argv[3]) # third argument for pdb filename with drude particle
#########################################################################################

## read/generate topology and force field
pdb = PDBFile(init_pdb)
pdb.topology.loadBondDefinitions(def_bonds_xml)
pdb.topology.createStandardBonds() # Create bonds based on the atom and residue names for all standard residue types. 
forcefield = ForceField(def_forces_xml)


## Add drude particle
## modeller: tools for editing molecular models
modeller = Modeller(pdb.topology, pdb.positions)
# add drude particles via forcefield xml file with positions 
#  which is the same position with previous atomtype in "def_forces_xml"
modeller.addExtraParticles(forcefield) 
#PDBFile.writeFile(modeller.topology, modeller.positions, open('init.pdb', 'w')) # check if correctly read


## construct OpenMM system representing a Topology with this force field.
system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=None, rigidWater=True)
# get memory access number for each forces described in "def_forces_xml"
nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]
customNonbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomNonbondedForce][0]
#custombond = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomBondForce][0]
# use basic pme method (return "4") for "NonbondedForce" part written in "def_forces_xml"
nbondedForce.setNonbondedMethod(NonbondedForce.PME) 
# use basic cut-off periodic method (return "2") of VdW LJ interaction for "customNonbondedForce" part written in "def_forces_xml"
customNonbondedForce.setNonbondedMethod(min(nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic)) 
print('nbMethod : ', customNonbondedForce.getNonbondedMethod()) # check if customNonbondedForce is assigned by VdW, not PME.
for i in range(system.getNumForces()):
  f = system.getForce(i)
  type(f)
  f.setForceGroup(i)
# you will see total 9 forces because the last CMMotionRemover is enabled by deafult with 8 forces are descrbied in "def_forces_xml"


## calculate total mass
totmass = 0.*dalton
for i in range(system.getNumParticles()):
	totmass += system.getParticleMass(i)
#system.addForce(AndersenThermostat(300*kelvin, 10/picosecond))


## platform settings
if 'gpu' in cpu_gpu_opencl:
	platform = Platform.getPlatformByName('CUDA')
	properties = {'CudaPrecision': 'mixed'}
elif 'cpu' in cpu_gpu_opencl:
	platform = Platform.getPlatformByName('CPU')
elif 'opencl' in cpu_gpu_opencl:
	platform = Platform.getPlatformByName('OpenCL')
	properties = {'OpenCLPrecision': 'mixed'}
else:
	raise ValueError(" wrong variable in cpu_gpu_opencl")


## integrator
integ_md = LangevinIntegrator(temperature, 1/picosecond, 0.001*picoseconds)
#integ_md = DrudeLangevinIntegrator(temperature, 1/picosecond, 1*kelvin, 1/picosecond, 0.001*picoseconds)
#barostat = MonteCarloBarostat(pressure,temperature,barofreq)
#system.addForce(barostat)
#barofreq = barostat.getFrequency()
#print(barofreq)


## re-positioning for drude particle from last point of pdb file
#pdb_new = PDBFile(last_pdb)
#modeller.topology.setUnitCellDimensions(pdb_new.topology.getUnitCellDimensions())
#modeller.positions = pdb_new.positions


## Simulation settings
## It ties together various objects used for running a simulation:
##   a Topology, System, Integrator, and Context. To use it, you provide the Topology, System, and Integrator,
##   and it creates the Context automatically.
if 'cpu' in cpu_gpu_opencl:
	simmd = Simulation(modeller.topology, system, integ_md, platform) # for CPU version
else:
	simmd = Simulation(modeller.topology, system, integ_md, platform, properties) # for GPU, OpenCL version
simmd.context.setPositions(modeller.positions)
platform =simmd.context.getPlatform()
platformname = platform.getName();
print(platformname)



#*************************************************************************************
#   This is section to create exclusions for TFSI nonbonded interactions, and update
#   Screened Drude interactions.  1-5 non-Coulomb interaction are accounted for
#   using CustomBondForce
#*************************************************************************************
print('Creating Exclusions for TFSI')

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
  else:
    for i in range(len(res._atoms)-1):
      for j in range(i+1,len(res._atoms)):
        (indi,indj) = (res._atoms[i].index, res._atoms[j].index)
      # here it doesn't matter if we already have this, since we pass the "True" flag
        nbondedForce.addException(indi,indj,0,1,0,True) # all intramolecular electrostatic interactions of TFSI are excluded
      # make sure we don't already exlude this customnonbond
        string1=str(indi)+"_"+str(indj)
        string2=str(indj)+"_"+str(indi)
        if string1 in flagexclusions and string2 in flagexclusions:
          continue
        else:
          customNonbondedForce.addExclusion(indi,indj) # all intramolecular VdW interactions are excluded

#************************************* end of TFSI section **************************
simmd.context.reinitialize()
simmd.context.setPositions(modeller.positions)

# load restart files
#simmd.loadCheckpoint('md_nvt.chk')

print('Starting calculation of cohesive energy.')

import mdtraj as md
# load trajectory
traj = md.load(cohesive_dcd,top=cohesive_top)
n_frames = len(traj.xyz)
for i_frame in range(n_frames):
  simmd.context.setPositions(traj.xyz[i_frame])
  state = simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
  if i_frame == n_frames - 1:
    # write initial pdb with drude oscillators
    position = state.getPositions()
    simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    PDBFile.writeFile(simmd.topology, position, open(strdir+'end_drudes.pdb', 'w'))
    for j in range(system.getNumForces()):
        f = system.getForce(j)
        print(type(f), str(simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()))
  print(str(state.getPotentialEnergy())) # print potential energy

print('Done!')

exit()
