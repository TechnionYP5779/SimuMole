from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
def scr(input_coor_name,simu_steps,temperature):
    pdb = app.PDBFile(input_coor_name)
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
         constraints=app.HBonds, rigidWater=True)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
        2.0*unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = mm.Platform.getPlatformByName('OpenCL')
    properties = {'OpenCLPrecision': 'mixed'}
    simulation = app.Simulation(pdb.topology, system, integrator, platform,
        properties)
    simulation.context.setPositions(pdb.positions)

    print('Minimizing...')
    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(temperature*unit.kelvin)
    print('Equilibrating...')
    simulation.step(100)

    simulation.reporters.append(app.DCDReporter('trajectory.dcd', 1000))
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
        speed=True, totalSteps=simu_steps, separator='\t'))

    print('Running Production...')
    simulation.step(simu_steps)
    print('Done!')