import pymol
from .simulation_main_script import create_movies_from_different_angles

temp = 'media/files/'  # path to temp folder


def pdb_and_dcd_match(pdb_file_name, dcd_file_name):  # todo: remove gui
    pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
    cmd = pymol.cmd
    try:
        cmd.reinitialize()
        cmd.load(temp + pdb_file_name)
        cmd.load_traj(temp + dcd_file_name)
        return True
    except:
        return False


def create_animations():
    # for now, only open pymol window: #todo: create the animations instead of only open PyMOL
    pymol.finish_launching(['pymol', '-q'])  # pymol: -q quiet launch, -c no gui, -e fullscreen
    cmd = pymol.cmd
    try:
        cmd.reinitialize()
        cmd.load(temp + "file_upload_pdb.pdb")
        cmd.load_traj(temp + "file_upload_dcd.dcd")
    except:  # mainly for "pymol.CmdException"
        update_simulation_status(
            'An error occurred while creating the simulation. Please try again later.')  # todo: change error message!
        return

    # create the animations:
    update_simulation_status('Creates the animations')
    create_movies_from_different_angles(cmd)  # create movies in media/movies folder

    # complete simulation:
    update_simulation_status('Done!')
    # self.cmd.quit() # todo: need to close PyMol window


def update_simulation_status(status):
    dir_path = temp
    simulation_status_path = dir_path + 'simulation_status.txt'
    with open(simulation_status_path, "w+") as f:
        f.write(status)
