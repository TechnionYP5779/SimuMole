from time import sleep
import pymol
import os


class Simulation:

    def __init__(self, first_pdb_id, second_pdb_id, x1, y1, z1, temperature):
        self.cmd = None

        self.first_pdb_id = first_pdb_id
        self.second_pdb_id = second_pdb_id
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.temperature = temperature

        # for debugging:
        self.first_pdb_id = '1GK7'
        self.second_pdb_id = '6CTH'

    def create_simulation(self):
        pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
        self.cmd = pymol.cmd

        # STEP 1: load input pdb
        self.save_pdb_by_id(self.first_pdb_id, 'first_pdb__' + str(self.first_pdb_id) + '.pdb')
        self.save_pdb_by_id(self.second_pdb_id, 'second_pdb__' + str(self.second_pdb_id) + '.pdb')

        # STEP 2: fix positions
        # (for now, only move the first protein. we also need to move the second)
        # todo: complete

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
        self.cmd.save('SimuMoleWeb/temp/' + name_of_file)
