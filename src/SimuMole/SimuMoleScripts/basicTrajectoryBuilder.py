from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout

dir_path = 'media/files/'
simulation_status_path = dir_path + 'simulation_status.txt'
simulation_status_during_run_path = dir_path + 'simulation_status_during_run.txt'
trajectory_path = dir_path + 'trajectory.dcd'


def scr(input_coor_name, simu_steps, temperature):
    pdb = app.PDBFile(input_coor_name)
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds,
                                     rigidWater=True)
    integrator = mm.LangevinIntegrator(temperature * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = mm.Platform.getPlatformByName('OpenCL')
    properties = {'OpenCLPrecision': 'mixed'}
    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)

    update_simulation_status('Minimizing energy...')
    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)

    update_simulation_status('Equilibrating...')
    simulation.step(100)

    simulation.reporters.append(app.DCDReporter(trajectory_path, 1000))
    simulation.reporters.append(
        app.StateDataReporter(simulation_status_during_run_path, 1000, progress=True, remainingTime=True,
                              totalSteps=simu_steps, separator=','))

    update_simulation_status('Running simulation...')
    simulation.step(simu_steps)
    print('Done!')

    update_simulation_status('Done!')


def update_simulation_status(status):
    with open(simulation_status_path, "w+") as f:
        f.write(status)


def scr_for_checks(input_coor_name):
    forcefield = app.ForceField('amber96.xml', 'tip3p.xml')
    pdb = app.PDBFile(input_coor_name)

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
                                     nonbondedCutoff=1.0 * unit.nanometers, constraints=app.HBonds, rigidWater=True,
                                     ewaldErrorTolerance=0.0005)
    integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds,
                                       2.0 * unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = mm.Platform.getPlatformByName('OpenCL')
    properties = {'OpenCLPrecision': 'mixed'}
    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)

    print('Minimizing...')
    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
    print('Equilibrating...')
    simulation.step(100)

    simulation.reporters.append(app.DCDReporter('media/files/' + 'very_good.dcd', 1000))
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True,
                                                      potentialEnergy=True, temperature=True, progress=True,
                                                      remainingTime=True,
                                                      speed=True, totalSteps=1000, separator='\t'))

    print('Running Production...')
    simulation.step(1000)
    print('Done!')
