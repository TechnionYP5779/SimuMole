from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit


def scr(input_coor_name, temperature, time_step_number, user_rand):
    dir_path = 'media/files/' + user_rand + '/'
    simulation_status_path = dir_path + 'simulation_status.txt'
    simulation_status_during_run_path = dir_path + 'simulation_status_during_run.txt'
    trajectory_path = dir_path
    traj_dcd = 'trajectory.dcd'

    pdb = app.PDBFile(input_coor_name)
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds,
                                     rigidWater=True)
    integrator = mm.LangevinIntegrator(temperature * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds)
    integrator.setConstraintTolerance(0.00001)

    platform = mm.Platform.getPlatformByName('CPU')

    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    update_simulation_status('Minimizing energy...', user_rand)
    simulation.minimizeEnergy()

    simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)

    update_simulation_status('Equilibrating...', user_rand)
    simulation.step(100)

    simulation.reporters.append(app.DCDReporter(trajectory_path + traj_dcd, 1000))
    simulation.reporters.append(
        app.StateDataReporter(simulation_status_during_run_path, 1000, progress=True, remainingTime=True,
                              totalSteps=time_step_number, separator=','))

    update_simulation_status('Running simulation...', user_rand)
    simulation.step(time_step_number)


def update_simulation_status(status, user_rand):
    dir_path = 'media/files/' + user_rand + '/'
    simulation_status_path = dir_path + 'simulation_status.txt'

    with open(simulation_status_path, "w+") as f:
        f.write(status)


def scr_for_checks(input_coor_name, user_rand, full_check=False):
    pdb = app.PDBFile(input_coor_name)
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

    # The following line causes an error (when we do not perform "fix pdb")
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds,
                                     rigidWater=True)

    if full_check:
        integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds)
        integrator.setConstraintTolerance(0.00001)

        platform = mm.Platform.getPlatformByName('CPU')

        simulation = app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)

        # Minimizing
        simulation.minimizeEnergy()

        simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

        # Equilibrating
        simulation.step(100)

        simulation.reporters.append(app.DCDReporter('media/files/' + user_rand + '/' + 'scr_for_checks.dcd', 1000))

        # Running Production
        simulation.step(1000)
