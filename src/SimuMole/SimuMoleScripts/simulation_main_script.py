from time import sleep
import pymol
import os

from .transformations import translate_pdb

temp = 'SimuMoleWeb/temp/'  # path to temp folder
pdb = '.pdb'  # pdb suffix


class Simulation:

    def __init__(self, first_pdb_id, second_pdb_id, x1, y1, z1, temperature):
        self.cmd = None

        self.first_pdb_id = first_pdb_id
        self.second_pdb_id = second_pdb_id
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.z1 = float(z1)
        # self.x2 = float(x2)
        # self.y2 = float(y2)
        # self.z2 = float(z2)
        self.temperature = float(temperature)

        # for debugging: # todo: delete when complete
        self.first_pdb_id, self.second_pdb_id = '1GK7', '6CTH'
        self.x1, self.y1, self.z1 = float(50), float(50), float(50)
        self.x2, self.y2, self.z2 = float(0), float(0), float(0)

    def create_simulation(self):
        pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
        self.cmd = pymol.cmd

        # STEP 1: load input pdb
        filename_1 = 'pdb_1__' + str(self.first_pdb_id)
        filename_2 = 'pdb_2__' + str(self.second_pdb_id)
        self.save_pdb_by_id(self.first_pdb_id, filename_1 + pdb)
        self.save_pdb_by_id(self.second_pdb_id, filename_2 + pdb)

        # STEP 2: fix positions
        filename_1_movement = filename_1 + '__movement'
        filename_2_movement = filename_2 + '__movement'
        translate_pdb(temp + filename_1 + pdb, temp + filename_1_movement + pdb, self.x1, self.y1, self.z1)
        translate_pdb(temp + filename_2 + pdb, temp + filename_2_movement + pdb, self.x2, self.y2, self.z2)

        # STEP 3: merge to single pdb file
        # todo: complete

        # STEP 4: use OpenMM
        # todo: complete

        # self.cmd.quit() # todo: need to close PyMol window

    def clear_simulation(self):  # todo: complete this! delete all temporary files
        # os.remove('path/to/files')
        return

    def save_pdb_by_id(self, pdb_id, name_of_file):
        self.cmd.reinitialize()
        sleep(0.5)

        self.cmd.load("https://files.rcsb.org/download/" + pdb_id + '.pdb')
        self.cmd.zoom()
        self.cmd.save(temp + name_of_file)