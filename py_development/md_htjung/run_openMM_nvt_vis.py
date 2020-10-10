from __future__ import print_function
import os
os.environ['OPENMM_CPU_THREADS'] = '6'
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime

temperature=323.15*kelvin
#pressure = Vec3(1.0,1.0,1.0)*atmosphere
#pressure = 1.0*atmosphere

pdb = PDBFile('nvt_init.pdb')
strdir = ''
cpu_gpu= 'GPU'

integ_md = DrudeLangevinIntegrator(temperature, 1/picosecond, 1*kelvin, 1/picosecond, 0.001*picoseconds)
integ_md.setMaxDrudeDistance(0.02)  # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)
print(" random seed = {}".format(integ_md.getRandomNumberSeed()))
integ_md.setRandomNumberSeed(1985)
print(" random seed = {}".format(integ_md.getRandomNumberSeed()))

pdb.topology.loadBondDefinitions('sapt_residues.xml')
pdb.topology.createStandardBonds();

modeller = Modeller(pdb.topology, pdb.positions)
forcefield = ForceField('sapt.xml')
modeller.addExtraParticles(forcefield)
#PDBFile.writeFile(modeller.topology, modeller.positions, open('init.pdb', 'w'))

system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=None, rigidWater=True, removeCMMotion=True)
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

#barofreq = 100
#barostat = MonteCarloBarostat(pressure,temperature,barofreq)
#system.addForce(barostat)

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
particleMap = {} # drude particle indice list
for i in range(drudeForce.getNumParticles()):
  particleMap[drudeForce.getParticleParameters(i)[0]] = i

# can't add duplicate ScreenedPairs, so store what we already have
flagexceptions = {} # all intramolecular 1-4 paris (including drudes) list
for i in range(nbondedForce.getNumExceptions()):
  (particle1, particle2, charge, sigma, epsilon) = nbondedForce.getExceptionParameters(i)
  string1=str(particle1)+"_"+str(particle2)
  string2=str(particle2)+"_"+str(particle1)
  flagexceptions[string1]=1
  flagexceptions[string2]=1
#print(len(flagexceptions))

# can't add duplicate customNonbonded exclusions, so store what we already have
flagexclusions = {} # all intramolecular paris (including drudes) list
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
# Actually, flagexclusion and flagexception dos not change by the for-loop above
# that is to just make sure

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
simmd.reporters.append(DCDReporter(strdir+'md_nvt.dcd',5,enforcePeriodicBox=False))

#*************************************************
# ChangYun created the DrudeDataReporter class, we need to pull this into this OpenMM install if we want to use it
#*************************************************
#simmd.reporters.append(DrudeDataReporter(strdir+'md_nvt.log', 1000, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, langevin=True, density=False,speed=True))
#simmd.reporters.append(DrudeDataReporter(strdir+'md_nvt_temp.log', 10000, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, langevin=True, drudeTemperature=True,density=False,speed=True))
#simmd.reporters.append(CheckpointReporter(strdir+'md_nvt.chk', 10000))
#simmd.reporters[1].report(simmd,state)
#simmd.reporters[2].report(simmd,state)

#for i in range(simmd.system.getNumForces()):
#    if type(simmd.system.getForce(i)) == MonteCarloBarostat:
#        simmd.system.removeForce(i)

print('Simulating...')
import numpy as np
def pairwise_row_diff_slicing_x(a):
  n = len(a)
  N = n*(n-1)//2
  idx = np.concatenate(( [0], np.arange(n-1,0,-1).cumsum() ))
  #print(idx)
  start, stop = idx[:-1], idx[1:]
  #print(start,stop)
  out = np.empty(N,dtype=a.dtype)
  for j,i in enumerate(range(n-1)):
    out[start[j]:stop[j]] = a[i+1:] - a[i,None]
  return out

def pairwise_row_diff_slicing_pbc_x(a,box,remove_idx):
  n = len(a)
  N = n*(n-1)//2
  idx = np.concatenate(( [0], np.arange(n-1,0,-1).cumsum() ))
  start, stop = idx[:-1], idx[1:]
  out = np.empty(N,dtype=a.dtype)
  for j,i in enumerate(range(n-1)):
    temp = a[i+1:] - a[i,None]
    out[start[j]:stop[j]] = temp - np.around(temp/box-0.5)*box
  out[remove_idx] = 0.
  return out

def read_pos_vel_for(state,remove_idx):
  position = state.getPositions(asNumpy=True)
  pos = position._value # nm
  velocity = state.getVelocities(asNumpy=True)
  vel = velocity.in_units_of(meter/second)._value # m/s
  forces = state.getForces(asNumpy=True)
  force = forces._value # kj/(nm*mole)
  #return np.delete(pos, remove_idx, 0), np.delete(vel, remove_idx, 0), np.delete(force, remove_idx, 0)
  return pos, vel, force

def kin_tensor(mass,vel):
  nvelx_e = vel[:,0]+0.0005*atomic_for[:,0]/mass[:] # 0.0005 = 0.5 * time step (ps)
  nvely_e = vel[:,1]+0.0005*atomic_for[:,1]/mass[:]
  kin_e = 0.5*np.sum(mass[:]*(nvelx_e[:]*nvely_e[:]))/1000 # KE_xy (kJ/mol)
  ## check total kinetic energy
  #for iatom in range(len(atomic_pos)):
  #    nvelx = vel[:,0]+0.0005*atomic_for[:,0]/mass[:]
  #    nvely = vel[:,1]+0.0005*atomic_for[:,1]/mass[:]
  #    nvelz = vel[:,2]+0.0005*atomic_for[:,2]/mass[:]
  #velx = 0.5*np.sum(mass*(nvelx**2))/1000
  #vely = 0.5*np.sum(mass*(nvely**2))/1000
  #velz = 0.5*np.sum(mass*(nvelz**2))/1000
  #kin_i= velx+vely+velz
  return kin_e

def virial_tensor(dist,pair_force):
  #return -np.sum(dist[:,0]*pair_force[:,1])/2. # xy
  return -np.sum(dist[:]*pair_force[:])/2. # xy

# generate mass array
mass_array = []
drude_particle = []
n_atoms = system.getNumParticles()
for i in range(n_atoms):
  imass = system.getParticleMass(i)
  if imass._value < 0.0: # drude particle
    drude_particle.append(i)
    continue
  mass_array.append(imass.in_units_of(kilogram/mole)._value)
mass_array = np.array(mass_array)

# make a indice list to exclude pairs for virial calculation based on exclusion list
n_atoms = system.getNumParticles()
n_pairs = n_atoms*(n_atoms-1)//2
idx_list = np.concatenate(( [0], np.arange(n_atoms-1,0,-1).cumsum() ))
virial_exclude = []
for i in range(nbondedForce.getNumExceptions()):
  (particle1, particle2, charge, sigma, epsilon) = nbondedForce.getExceptionParameters(i)
  minv = np.min((particle1,particle2))
  maxv = np.max((particle1,particle2))
  virial_exclude.append(int(idx_list[minv]+maxv-minv-1))
virial_exclude = np.array(virial_exclude)
#print(len(virial_exclude))

for i in range(1,11):
  ### only for Verlet velocity and Leapfrog integrator ############
  ## read velocity just before (dt-1/2)
  #simmd.step(4)
  #state = simmd.context.getState(getVelocities=True)
  #velocity_prev = state.getVelocities(asNumpy=True)
  #atomic_vel_prev = velocity_prev.in_units_of(meter/second)._value # m/s
  #atomic_vel_prev = np.delete(atomic_vel_prev, drude_particle, 0)
  ##################################################################
  simmd.step(5)
  
  # read velocity at the time (dt+1/2)
  state = simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,enforcePeriodicBox=True)
  atomic_pos, atomic_vel, atomic_for = read_pos_vel_for(state,drude_particle)
  
  # kinetic energy term
  kin_term = kin_tensor(mass_array,atomic_vel)
  print(" kin. energy (xy) by mannual= {} kJ/mol".format(kin_term))
  #print(" kin. energy by OpenMM = {}".format(str(state.getKineticEnergy()))) 

  ## check center of mass velocity
  #reduce_atomic_vel = np.empty(np.shape(atomic_vel))
  #n_atoms = len(mass_array)
  #for i_atom in range(n_atoms):
  #  reduce_atomic_vel[i_atom] = mass_array[i_atom] * atomic_vel[i_atom]
  #com_vel = np.sum(reduce_atomic_vel, axis=0)/np.sum(mass_array)
  #com_kin_e = 0.5*np.sum(mass_array)*np.sum(com_vel*com_vel)
  #if com_kin_e*100/total_kin_e > 1.0:
  #  print(" COM kinetic energy = {} kJ/mol (should be removed!)".format(com_kin_e))
  #  print(total_kin_e - com_kin_e)

  # virial term
  box_vec = state.getPeriodicBoxVectors()
  box_3dim = box_vec._value[0][0],  box_vec._value[1][1], box_vec._value[2][2]
  print("{} {} max".format(np.amax(atomic_for[:,1]),np.argmax(atomic_for[:,1])))
  pair_dfor = pairwise_row_diff_slicing_x(atomic_for[:,1]) # difference between any two vectors in atomic_pos without duplicates
  #print("{} {} max".format(np.amax(pair_dfor),np.argmax(pair_dfor)))
   
  pair_dis_short = pairwise_row_diff_slicing_pbc_x(atomic_pos[:,0],box_3dim[0],virial_exclude) 
  vir_term = virial_tensor(pair_dis_short,pair_dfor)
  #print(" virial energy (xy) = {} kJ/mol".format(vir_term))

  # pressure tensor
  ####box_vol = state.getPeriodicBoxVolume() # nm^3
  ####box_vol = box_vol._value
  ####pair_dfor0 = pairwise_row_diff_slicing_x(atomic_for[:,0])
  ####pair_dfor1 = pairwise_row_diff_slicing_x(atomic_for[:,1])
  ####pair_dfor2 = pairwise_row_diff_slicing_x(atomic_for[:,2])
  ####pair_dis_short0 = pairwise_row_diff_slicing_pbc_x(atomic_pos[:,0],box_3dim[0],virial_exclude) 
  ####pair_dis_short1 = pairwise_row_diff_slicing_pbc_x(atomic_pos[:,1],box_3dim[1],virial_exclude) 
  ####pair_dis_short2 = pairwise_row_diff_slicing_pbc_x(atomic_pos[:,2],box_3dim[2],virial_exclude) 
  ####vir_term0 = virial_tensor(pair_dis_short0,pair_dfor0)
  ####print(np.argmax(vir_term0))
  ####print("{} {} {} max".format(np.amax(pair_dfor0),pair_dis_short0[np.argmax(pair_dfor0)],pair_dfor0[np.argmax(pair_dfor0)]))
  ####print("{} {} {} min".format(np.amin(pair_dfor0),pair_dis_short0[np.argmin(pair_dfor0)],pair_dfor0[np.argmin(pair_dfor0)]))
  ####vir_term1 = virial_tensor(pair_dis_short1,pair_dfor1)
  ####vir_term2 = virial_tensor(pair_dis_short2,pair_dfor2)
  ####vir_term_iso = (vir_term0+vir_term1+vir_term2)/3.
  ####pressure_tensor = 2.*(kin_term - vir_term)*100/(box_vol*6.0221417930) # unit: bar
  ####pressure_tensor_iso = 2.*(state.getKineticEnergy()._value - vir_term_iso)*100/(box_vol*6.0221417930) # unit: bar
  #####print("est. pressure (iso) = {} bar".format((pressure_tensor[0]+pressure_tensor[3]+pressure_tensor[5])/3.))
  ####print(" p= {}".format(pressure_tensor))
  ####print(" P(iso)= {}".format(pressure_tensor_iso))
  #print(i,strftime("%Y-%m-%d %H:%M:%S", gmtime()))
  #print(i,datetime.now())
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
PDBFile.writeFile(simmd.topology, position, open(strdir+'md_nvt.pdb', 'w'))

np.savetxt('pres.out',pressure_tensor_out)
np.save('pres.out',pressure_tensor_out)

print('Done!')

exit()
