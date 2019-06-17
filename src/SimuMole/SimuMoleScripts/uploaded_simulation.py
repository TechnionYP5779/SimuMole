import pymol2
from .simulation_main_script import create_movies_from_different_angles
from .basicTrajectoryBuilder import update_simulation_status
from time import sleep

temp = 'media/files/'  # path to temp folder


def pdb_and_dcd_match(pdb_file_name, dcd_file_name, user_rand):
    # pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
    # cmd = pymol.cmd
    p1 = pymol2.PyMOL()
    p1.start()
    cmd = p1.cmd

    try:
        cmd.reinitialize()
        sleep(0.5)
        cmd.load(temp + user_rand + '/' + pdb_file_name)
        cmd.load_traj(temp + user_rand + '/' + dcd_file_name)
        return True
    except Exception as e:
        print(str(e))
        return False


def create_animations(user_rand):
    try:
        # pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
        # cmd = pymol.cmd
        p1 = pymol2.PyMOL()
        p1.start()
        cmd = p1.cmd

        # load input
        cmd.reinitialize()
        sleep(0.5)
        cmd.load(temp + user_rand + "/file_upload_pdb.pdb")
        cmd.load_traj(temp + user_rand + "/file_upload_dcd.dcd")

        # create the animations:
        update_simulation_status('Creates the animations', user_rand)
        create_movies_from_different_angles(cmd, user_rand)  # create movies in media/movies folder

        # complete simulation:
        update_simulation_status('Done!', user_rand)

    except Exception as e:  # mainly for "pymol.CmdException"
        print(str(e))
        update_simulation_status('An error occurred while creating the animations. Please try again later.', user_rand)
        return
