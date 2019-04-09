
from time import sleep
import pymol
pymol.finish_launching(['pymol', '-qc'])  # pymol: quiet and no GUI

for index in [1, 2, 3]:
    #  If using a single script to produce multiple figures,
    #  then the reinitialize() command is necessary between parts of the script
    pymol.cmd.reinitialize()
    # Desired pymol commands here to produce and save figures
    sleep(0.5)  # (in seconds)
