from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

'''
This script utilizes the pdbfixer package in order to automatically fix pdb files for OpenMM.
If the pdbfixer import doesn't work, install the package with the following command through the Pycharm terminal:

conda install -c omnia -c conda-forge pdbfixer


The script itself takes a pdb file name as an argument and rewrites it.
'''


def fix_pdb(pdb_file):
    fixer = PDBFixer(filename=pdb_file)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdb_file, 'w'))
