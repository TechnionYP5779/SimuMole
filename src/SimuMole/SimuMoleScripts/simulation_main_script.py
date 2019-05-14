from time import sleep
import pymol
import os

from .fix_pdb import fix_pdb
from .basicTrajectoryBuilder import scr
from .transformations import translate_pdb
from .OpenMM_scriptBuilder import create_openmm_script, openMMbuilder

temp = 'media/files/'  # path to temp folder
pdb = '.pdb'  # pdb suffix


class Simulation:

    def __init__(self, num_of_proteins, first_pdb_type, first_pdb_id, second_pdb_type, second_pdb_id,
                 x1, y1, z1, x2, y2, z2, degXY_1, degYZ_1, degXY_2, degYZ_2, temperature, production_steps):

        self.cmd = None

        self.num_of_proteins = num_of_proteins

        self.first_pdb_type = first_pdb_type
        self.first_pdb_id = first_pdb_id

        self.second_pdb_type = second_pdb_type
        self.second_pdb_id = second_pdb_id
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.z1 = float(z1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.z2 = float(z2)
        self.degXY_1 = float(degXY_1)
        self.degYZ_1 = float(degYZ_1)
        self.degXY_2 = float(degXY_2)
        self.degYZ_2 = float(degYZ_2)
        self.temperature = float(temperature)
        self.production_steps = int(production_steps)

    #   # for debugging: # todo: delete when complete with debugging
    #   self.first_pdb_id, self.second_pdb_id = '1GK7', '6CTH'
    #   self.x1, self.y1, self.z1 = float(50), float(50), float(50)
    #   self.x2, self.y2, self.z2 = float(0), float(0), float(0)

    def create_simulation(self):
        pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
        self.cmd = pymol.cmd

        filename_1 = '_1_'
        filename_2 = '_2_'

        if self.num_of_proteins == '2':
            # STEP 1: load input pdb
            if self.first_pdb_type == 'by_id':
                self.save_pdb_by_id(self.first_pdb_id, filename_1 + pdb)
            if self.second_pdb_type == 'by_id':
                self.save_pdb_by_id(self.second_pdb_id, filename_2 + pdb)

            # STEP 2: fix positions
            filename_1_movement = filename_1 + '__movement'
            filename_2_movement = filename_2 + '__movement'
            translate_pdb(temp + filename_1 + pdb, temp + filename_1_movement + pdb, self.x1, self.y1, self.z1,
                          self.degXY_1, self.degYZ_1)
            translate_pdb(temp + filename_2 + pdb, temp + filename_2_movement + pdb, self.x2, self.y2, self.z2,
                          self.degXY_2, self.degYZ_2)

            # STEP 2.5: fix pdb
            fix_pdb(temp + filename_1_movement + pdb)
            fix_pdb(temp + filename_2_movement + pdb)

            # STEP 3: merge to single pdb file
            self.save_pdbs_in_one_pdb(filename_1_movement, filename_2_movement)

            # STEP 3.5: fix pdb
            fix_pdb(temp + "both__" + filename_1_movement + '_' + filename_2_movement + pdb)

            # STEP 4: use OpenMM # todo: complete
            input_coor_name = temp + "both__" + filename_1_movement + '_' + filename_2_movement + pdb

           # openMMbuilder('media/files/', input_coor=input_coor_name, state_dataT=True, pdbT=True, dcdT=False,
           #       report_interval=1000, equilibration_steps=100, production_steps=1000, minimize=True,
           #       max_minimize_steps=3, temperature=self.temperature, platform='OpenCL')

            scr(input_coor_name,self.production_steps,self.temperature)


        else:
            # STEP 1: load input pdb
            if self.first_pdb_type == 'by_id':
                self.save_pdb_by_id(self.first_pdb_id, filename_1 + pdb)

            # STEP 2: fix positions
            filename_1_movement = filename_1 + '__movement'
            translate_pdb(temp + filename_1 + pdb, temp + filename_1_movement + pdb, self.x1, self.y1, self.z1)

            # STEP 2.5: fix pdb
            fix_pdb(temp + filename_1_movement + pdb)

            # STEP 3: use OpenMM # todo: complete
            input_coor_name = temp + filename_1_movement + pdb
            # openMMbuilder('media/files/', input_coor=input_coor_name, state_dataT=True, pdbT=True, dcdT=False,
            #              report_interval=1000, equilibration_steps=100, production_steps=1000, minimize=True,
            #              max_minimize_steps=1, temperature=self.temperature, platform='OpenCL')


            scr(input_coor_name, self.production_steps, self.temperature) 




        self.cmd.reinitialize()
        self.cmd.load(input_coor_name)
        self.cmd.load(temp + 'trajectory.dcd')
        # self.cmd.quit() # todo: need to close PyMol window

    def clear_simulation(self):  # todo: complete this! delete all temporary files
        # os.remove('path/to/files')
        return

    def save_pdb_by_id(self, pdb_id, name_of_file):
        self.cmd.reinitialize()
        sleep(0.5)

        self.cmd.load("https://files.rcsb.org/download/" + pdb_id + pdb)
        self.cmd.zoom()
        self.cmd.save(temp + name_of_file)

    def save_pdbs_in_one_pdb(self, filename_1, filename_2):
        self.cmd.reinitialize()
        sleep(0.5)

        self.cmd.load(temp + filename_1 + pdb)
        self.cmd.load(temp + filename_2 + pdb)
        self.cmd.zoom()
        self.cmd.save(temp + "both__" + filename_1 + '_' + filename_2 + pdb)

    @staticmethod
    def merge_pdbs_by_copy(filename_1, filename_2, do_fix=True):
        merged_pdb = temp + "both__" + filename_1 + "_" + filename_2 + pdb
        merged_file = open(merged_pdb, 'w')
        file1 = open(temp + filename_1 + pdb)
        file2 = open(temp + filename_2 + pdb)
        for line in file1:
            if not (line.startswith('MASTER') or line.startswith('END')):
                merged_file.write(line)
        for line in file2:
            if not (line.startswith('HEADER')):
                merged_file.write(line)
        file1.close()
        file2.close()
        merged_file.close()
        if do_fix:
            fix_pdb(merged_pdb)
