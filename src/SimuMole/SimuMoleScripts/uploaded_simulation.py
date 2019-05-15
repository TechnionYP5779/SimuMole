import pymol
import os
temp = 'media/files/'  # path to temp folder

class Uploaded_Simulation:
    def __init__(self, pdb_file_name, dcd_file_name):
        self.cmd = None
        self.pdb_file_name = pdb_file_name
        self.dcd_file_name = dcd_file_name

    def run_simulation(self):
        pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
        self.cmd = pymol.cmd
        self.cmd.load(temp + self.pdb_file_name)
        self.cmd.load(temp + self.dcd_file_name)