#import pymol
import numpy as np


COORDINATES_START = 31
COORDINATES_LEN = 7


'''
Gathers the coordinates of all atoms in a pdb file.
Arguments:
pdb - the name of a pdb file.
Returns:
A list of all coordinates, each in a numpy array of size 3.
'''
def get_atoms(pdb):
    file = open(pdb)
    vecs = []
    for line in file:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        words = list(filter(lambda s: s != '', line[COORDINATES_START:].split(' ')))
        vecs.append(np.array([float(words[0]), float(words[1]), float(words[2])]))
    file.close()
    return vecs


'''
Prepares a number to string form for the pdb file.
Due to the limitations of the format, each coordinate can have no more than 6 characters not including the minus.
Thus this function takes a number and changes the string to make it easier to insert into the pdb file.
Aguments:
n - the number to be prepared.
Returns:
A string of the prepared number. The string will always be 7 characters wide.
'''
def prepare_num(n):
    result = str(n)
    if n >= 0:
        result = ' ' + result
    result = result.ljust(COORDINATES_LEN, ' ')
    return result[:COORDINATES_LEN]


'''
Takes the coordinates of an atom and outputs a string that can be inserted into a pdb file.
'''
def prepare_coords(coords):
    return prepare_num(coords[0]) + ' ' + prepare_num(coords[1]) + ' ' + prepare_num(coords[2])


'''
Takes an existing pdb file and replaces the coordinates of its atoms.
If there are not enough coordinates to place into the file, the remaining coordinates will be replaced with zeros.
Arguments:
old_pdb - the name of the old pdb file.
new_pdb - the name of the new file that will be written by the function.
vecs - the new coordinates.
'''
def change_pdb(old_pdb, new_pdb, vecs):
    if old_pdb == new_pdb:
        return
    old_file = open(old_pdb)
    new_file = open(new_pdb, 'w')
    cnt = 0
    for line in old_file:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            new_file.write(line)
        else:
            new_coord = np.array([0, 0, 0])
            if cnt < len(vecs):
                new_coord = vecs[cnt]
                cnt += 1
            new_line = line[:31] + prepare_coords(new_coord) + line[54:]
            new_file.write(new_line)
    old_file.close()
    new_file.close()


'''
Applies translation transformation to the atoms coordinates, moving every atom a specified distance in space.
Arguments:
x - the value added to each x coordinate.
y - the value added to each y coordinate.
z - the value added to each z coordinate.
vecs - list of coordinates, this list will be changed by this function.
'''
def translate_vecs(x, y, z, vecs):
    for v in vecs:
        v[0] += x
        v[1] += y
        v[2] += z


'''
Applies translation transformation to a pdb file, resulting in a new pdb file with different atom coordinates.
Arguments:
old_pdb - the name of the old pdb file.
new_pdb - the name of the new file that will be written by the function.
x - the amount each atom will be moved in the X axis.
y - the amount each atom will be moved in the Y axis.
z - the amount each atom will be moved in the Z axis.
'''
def translate_pdb(old_pdb, new_pdb, x, y, z):
    vecs = get_atoms(old_pdb)
    translate_vecs(x, y, z, vecs)
    change_pdb(old_pdb, new_pdb, vecs)
