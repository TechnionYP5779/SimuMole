
from time import sleep
import pymol_plugin_dynamics
pymol_plugin_dynamics.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen

cmd = pymol_plugin_dynamics.cmd

for index in [1, 2, 3]:
    #  If using a single script to produce multiple figures,
    #  then the reinitialize() command is necessary between parts of the script
    cmd.reinitialize()
    # Desired pymol commands here to produce and save figures
    sleep(0.5)  # (in seconds)

    # cmd.load("BSAP_one_chain.pdb")
    # cmd.load("model1asn.pdb")


def pdb_operations(local_path_a=None, local_path_b=None, id_a=None, id_b=None):
    if local_path_b is None:
        if local_path_a is not None:
            cmd.load(local_path_a)
    else:
        cmd.load(local_path_a)
        cmd.load(local_path_b)

    if id_b is None:
        if id_a is not None:
            cmd.load("https://files.rcsb.org/download/" + id_a)
    else:
        cmd.load("https://files.rcsb.org/download/" + id_a)
        cmd.load("https://files.rcsb.org/download/" + id_b)

    cmd.zoom()  # defaults to zoom = all
    cmd.save("all.pdb")


#  pdb_operations("BSAP_one_chain.pdb", "model1asn.pdb")

