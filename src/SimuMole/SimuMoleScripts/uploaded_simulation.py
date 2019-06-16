import pymol
from .simulation_main_script import create_movies_from_different_angles
from .basicTrajectoryBuilder import update_simulation_status
from time import sleep

temp = 'media/files/'  # path to temp folder


def pdb_and_dcd_match(pdb_file_name, dcd_file_name):
    pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
    cmd = pymol.cmd
    try:
        cmd.reinitialize()
        sleep(0.5)
        cmd.load(temp + pdb_file_name)
        cmd.load_traj(temp + dcd_file_name)
        return True
    except Exception as e:
        print(str(e))
        return False


def create_animations():
    try:
        pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
        cmd = pymol.cmd

        # load input
        cmd.reinitialize()
        sleep(0.5)
        cmd.load(temp + "file_upload_pdb.pdb")
        cmd.load_traj(temp + "file_upload_dcd.dcd")

        # create the animations:
        update_simulation_status('Creates the animations')
        create_movies_from_different_angles(cmd)  # create movies in media/movies folder

        # complete simulation:
        update_simulation_status('Done!')

    except Exception as e:  # mainly for "pymol.CmdException"
        print(str(e))
        update_simulation_status('An error occurred while creating the animations. Please try again later.')
        return
