import pymol

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
    cmd.reinitialize()
    cmd.load(temp + "file_upload_pdb.pdb")
    cmd.load(temp + "file_upload_dcd.dcd")
