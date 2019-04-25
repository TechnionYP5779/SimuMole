from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

pdb = app.PDBFile('both__pdb_1__6DY4__movement_pdb_2__5ZKH__movement.pdb')
forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, 
     constraints=app.HBonds, rigidWater=True)
integrator = mm.LangevinIntegrator(20.0*unit.kelvin, 1.0/unit.picoseconds, 
    2.0*unit.femtoseconds)
integrator.setConstraintTolerance(1e-05)

platform = mm.Platform.getPlatformByName('OpenCL')
properties = {'OpenCLPrecision': 'mixed'}
simulation = app.Simulation(pdb.topology, system, integrator, platform, 
    properties)
simulation.context.setPositions(pdb.positions)

print('Minimizing...')
simulation.minimizeEnergy(maxIterations=3)

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
print('Equilibrating...')
simulation.step(100)

simulation.reporters.append(app.PDBReporter('trajectory.pdb', 1000))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, 
time=False, potentialEnergy=True, kineticEnergy=False, totalEnergy=False, 
temperature=True, volume=False, density=False, progress=True, remainingTime=True, speed=True, totalSteps=1000, separator='\t'))

print('Running Production...')
simulation.step(1000)
print('Done!')
