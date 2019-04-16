import io

"""
The script below generate the openmm code.
The input is the parameters listed below - each parameter has its own default value. I used the default measures too.
There are 2 options:
	A - create a file named 'openmm.py' with the generated code
	B - run the generated code
You can choose between the two options in the beginning of the function
"""


def openMMbuilder(output_path='', input_coor='input.pdb', force_field='amber99sbildn', water_model='tip3p',
                  platform='CUDA',
                  precision='mixed', device_index=-1,
                  OpenCL_platform_index=-1, nonbonded_method='PME', ewald_error_tolerance=0.0005, constraints='HBonds',
                  constraint_error_tol=0.00001,
                  rigid_water=True, nonbonded_cutoff=1.0, random_init_vels=True, generation_temp=300,
                  integrator='Langevin', time_step=2.0, error_tolerance=0.0001, collision_rate=1.0, temperature=300,
                  barostat='None', pressure=1, barostat_interval=25,
                  state_dataT=True, dcdT=True, pdbT=False, report_interval=1000, equilibration_steps=100,
                  production_steps=1000, minimize=True,
                  max_minimize_steps=-1,
                  state_data_options=[True, False, True, True, True, False, False, True, False, False]):
    # for option A put True in optionAB, for option B put False
    optionAB = True

    code = ''

    code += ('from __future__ import print_function\n')
    code += ('from simtk.openmm import app\n')
    code += ('import simtk.openmm as mm\n')
    code += ('from simtk import unit\n')
    code += ('from sys import stdout\n')

    code += ('\n')

    code += ('pdb = app.PDBFile(\'' + input_coor + '\')\n')
    code += ('forcefield = app.ForceField(\'' + force_field + '.xml\', \'' + water_model + '.xml\')\n')

    code += ('\n')

    if constraints != 'None':
        constraints_tmp = 'constraints=app.' + constraints + ', '
    else:
        constraints_tmp = 'constraints=None, '

    if nonbonded_method == 'PME' or nonbonded_method == 'Ewald':
        code += ('system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.' + nonbonded_method + ', \n    ' +
                 'nonbondedCutoff=' + str(
                    nonbonded_cutoff) + '*unit.nanometers, ' + constraints_tmp + 'rigidWater=' + str(
                    rigid_water) + ', \n    ' +
                 'ewaldErrorTolerance=' + str(ewald_error_tolerance) + ')\n')
    elif nonbonded_method == 'NoCutoff':
        code += (
                'system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.' + nonbonded_method + ', \n     ' +
                constraints_tmp + 'rigidWater=' + str(rigid_water) + ')\n')
    else:
        code += ('system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.' + nonbonded_method + ', \n    ' +
                 'nonbondedCutoff=' + str(
                    nonbonded_cutoff) + '*unit.nanometers, ' + constraints_tmp + 'rigidWater=' + str(
                    rigid_water) + ')\n')

    if integrator == 'Langevin':
        code += ('integrator = mm.LangevinIntegrator(' + str(temperature) + '*unit.kelvin, ' + str(
            collision_rate) + '/unit.picoseconds, \n    ' +
                 str(time_step) + '*unit.femtoseconds)\n')
    elif integrator == 'Verlet':
        code += ('integrator = mm.VerletIntegrator(' + str(time_step) + '*unit.femtoseconds)\n')
    elif integrator == 'Brownian':
        code += ('integrator = mm.BrownianIntegrator(' + str(temperature) + '*unit.kelvin, ' + str(
            collision_rate) + '/unit.picoseconds, \n    ' +
                 str(time_step) + '*unit.femtoseconds)\n')
    elif integrator == 'VariableLangevin':
        code += ('mm.VariableLangevinIntegrator(' + str(error_tolerance) + ')\n')
    elif integrator == 'VariableVerlet':
        code += ('integrator = mm.VariableVerletIntegrator(' + str(error_tolerance) + ')\n')

    if constraints != 'None':
        code += ('integrator.setConstraintTolerance(' + str(constraint_error_tol) + ')\n')

    if (barostat != 'None'):
        code += ('system.addForce(mm.MonteCarloBarostat(' + str(pressure) + '*unit.atmospheres, ' + str(
            temperature) + '*unit.kelvin, ' + str(barostat_interval) + '))\n')

    code += ('\n')

    if platform == 'CUDA' or platform == 'OpenCL':
        code += ('platform = mm.Platform.getPlatformByName(\'' + platform + '\')\n')
        if platform == 'CUDA':
            if device_index == -1:
                code += ('properties = {\'CudaPrecision\': \'' + precision + '\'}\n')
            else:
                code += ('properties = {\'CudaPrecision\': \'' + precision + '\', \'CudaDeviceIndex\': \'' + str(
                    device_index) + '\'}\n')

        else:
            if device_index == -1 and OpenCL_platform_index == -1:
                code += ('properties = {\'OpenCLPrecision\': \'' + precision + '\'}\n')
            else:
                if device_index == -1:
                    code += (
                            'properties = {\'OpenCLPrecision\': \'' + precision + '\', \'OpenCLPlatformIndex\': \'' + str(
                        OpenCL_platform_index) + '\'}\n')
                elif OpenCL_platform_index == -1:
                    code += (
                            'properties = {\'OpenCLPrecision\': \'' + precision + '\', \'OpenCLDeviceIndex\': \'' + str(
                        device_index) + '\'}\n')
                else:
                    code += (
                            'properties = {\'OpenCLPrecision\': \'' + precision + '\', \'OpenCLPlatformIndex\': \'' + str(
                        OpenCL_platform_index) + '\', \n              '
                            + '\'OpenCLDeviceIndex\': \'' + str(device_index) + '\'}\n')
        code += ('simulation = app.Simulation(pdb.topology, system, integrator, platform, \n    ' +
                 'properties)\n')

    else:
        code += ('platform = mm.Platform.getPlatformByName(\'' + platform + '\')\n')
        code += ('simulation = app.Simulation(pdb.topology, system, integrator, platform)\n')
    code += ('simulation.context.setPositions(pdb.positions)\n')

    code += ('\n')

    if minimize:
        code += ('print(\'Minimizing...\')\n')
        if max_minimize_steps > -1:
            code += ('simulation.minimizeEnergy(maxIterations=' + str(max_minimize_steps) + ')\n')
        else:
            code += ('simulation.minimizeEnergy()\n')

    code += ('\n')

    if random_init_vels:
        code += ('simulation.context.setVelocitiesToTemperature(' + str(generation_temp) + '*unit.kelvin)\n')
    code += ('print(\'Equilibrating...\')\n')
    code += ('simulation.step(' + str(equilibration_steps) + ')\n')

    code += ('\n')

    if dcdT:
        code += ('simulation.reporters.append(app.DCDReporter(\'trajectory.dcd\', ' + str(report_interval) + '))\n')

    if pdbT:
        code += ('simulation.reporters.append(app.PDBReporter(\'trajectory.pdb\', ' + str(report_interval) + '))\n')

    if state_dataT:
        code += ('simulation.reporters.append(app.StateDataReporter(stdout, ' + str(report_interval) + ', step=' + str(
            str(state_data_options[0])) + ', \n' +
                 'time=' + str(state_data_options[1]) + ', potentialEnergy=' + str(
                    state_data_options[4]) + ', kineticEnergy=' + str(state_data_options[5]) +
                 ', totalEnergy=' + str(state_data_options[6]) + ', \n' + 'temperature=' + str(state_data_options[7]) +
                 ', volume=' + str(state_data_options[8]) + ', density=' + str(state_data_options[9]) +
                 ', progress=' + str(state_data_options[3]) + ', remainingTime=' + str(state_data_options[3])
                 + ', speed=' + str(state_data_options[2]) + ', totalSteps=' + str(
                    production_steps) + ', separator=\'\\t\'))\n')

    code += ('\n')

    code += ('print(\'Running Production...\')\n')
    code += ('simulation.step(' + str(production_steps) + ')\n')
    code += ('print(\'Done!\')\n')

    # option 1: create a file name 'openmm.py' with the generated script
    if optionAB == True:
        f = open(output_path + 'openmm.py', 'w')
        f.write(code)

        # option 2: execute the code generated
        if optionAB == False:
            exec(code)


def create_openmm_script(first_pdb_id, second_pdb_id):
    input_coor_name = '' + "both__" + first_pdb_id + '_' + second_pdb_id + '.pdb'

    openMMbuilder('../SimuMoleWeb/temp/', input_coor=input_coor_name, state_dataT=True, pdbT=True, dcdT=False,
                  report_interval=11, equilibration_steps=10001, production_steps=5000, minimize=True,
                  max_minimize_steps=3)

# first_pdb_id, second_pdb_id = '1GK7', '6CTH'
# create_openmm_script(first_pdb_id, second_pdb_id)
