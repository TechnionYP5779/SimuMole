import numpy as np
import math

COORDINATES_START = 31
COORDINATES_LEN = 7


def get_atoms(pdb):
    """
    Gathers the coordinates of all atoms in a pdb file.
    :param pdb: the name of a pdb file.
    :return: A list of all coordinates, each in a numpy array of size 3.
    """
    file = open(pdb)
    vecs = []
    for line in file:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        words = list(filter(lambda s: s != '', line[COORDINATES_START:].split(' ')))
        vecs.append(np.array([float(words[0]), float(words[1]), float(words[2])]))
    file.close()
    return vecs


def get_atoms_string(str):
    """
    Gathers the coordinates of all atoms in string.
    :param str: the string containing the atoms, each atom separated by new line
    :return: A list of all coordinates, each in a numpy array of size 3.
    """
    vecs = []
    for line in str.splitlines():
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        words = list(filter(lambda s: s != '', line[COORDINATES_START:].split(' ')))
        vecs.append(np.array([float(words[0]), float(words[1]), float(words[2])]))
    return vecs


def prepare_num(n):
    """
    Prepares a number to string form for the pdb file.
    Due to the limitations of the format, each coordinate can have no more than 6 characters not including the minus.
    Thus this function takes a number and changes the string to make it easier to insert into the pdb file.
    :param n: the number to be prepared.
    :return: A string of the prepared number. The string will always be 7 characters wide.
    """
    result = str(n)
    if n >= 0:
        result = ' ' + result
    result = result.ljust(COORDINATES_LEN, ' ')
    return result[:COORDINATES_LEN]


def prepare_coords(coords):
    """
    Takes the coordinates of an atom and outputs a string that can be inserted into a pdb file.
    """
    return prepare_num(coords[0]) + ' ' + prepare_num(coords[1]) + ' ' + prepare_num(coords[2])


def change_pdb(old_pdb, new_pdb, vecs):
    """
    Takes an existing pdb file and replaces the coordinates of its atoms.
    If there are not enough coordinates to place into the file, the remaining coordinates will be replaced with zeros.
    :param old_pdb: the name of the old pdb file.
    :param new_pdb: the name of the new file that will be written by the function.
    :param vecs: the new coordinates.
    :return:
    """
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


def translate_vecs(x, y, z, vecs):
    """
    Applies translation transformation to the atoms coordinates, moving every atom a specified distance in space.
    :param x: the value added to each x coordinate.
    :param y:  the value added to each y coordinate.
    :param z: the value added to each z coordinate.
    :param vecs: list of coordinates, this list will be changed by this function.
    :return:
    """
    for v in vecs:
        v[0] += x
        v[1] += y
        v[2] += z


def translate_pdb(old_pdb, new_pdb, x, y, z, degXY, degYZ):
    """
    Applies translation transformation to a pdb file, resulting in a new pdb file with different atom coordinates.
    :param old_pdb: the name of the old pdb file.
    :param new_pdb: the name of the new file that will be written by the function.
    :param x: the amount each atom will be moved in the X axis.
    :param y: the amount each atom will be moved in the Y axis.
    :param z: the amount each atom will be moved in the Z axis.
    :param degXY: the rotation rate in degrees in XY level
    :param degYZ: the rotation rate in degrees in YZ level
    """
    vecs = get_atoms(old_pdb)
    translate_vecs(x, y, z, vecs)
    rotate_molecular(x, y, z, degXY, degYZ, vecs)
    change_pdb(old_pdb, new_pdb, vecs)


def rotate_molecular(x, y, z, degXY, degYZ, vecs):
    """
    Rotate each point in vecs in the degrees according to degXZ and degYZ around the center (x, y, z)
    :param x: the x axis of the center
    :param y: the y axis of the center
    :param z: the z axis of the center
    :param degXY: the rotation rate in degrees in XY level
    :param degYZ: the rotation rate in degrees in YZ level
    :param vecs: the vector of the points to be rotated
    Note: In the future we may calculate the center instead of getting it as parameters
    """
    degXY_Rad = math.radians(degXY)
    degYZ_Rad = math.radians(degYZ)
    center = x, y, z
    oX, oY, oZ = center

    for v in vecs:
        # rotation in XY level
        pX, pY = v[0], v[1]
        v[0] = oX + math.cos(degXY_Rad) * (pX - oX) - math.sin(degXY_Rad) * (pY - oY)
        v[1] = oY + math.sin(degXY_Rad) * (pX - oX) + math.cos(degXY_Rad) * (pY - oY)

        # rotation in YZ level
        pY, pZ = v[1], v[2]
        v[1] = oY + math.cos(degYZ_Rad) * (pY - oY) - math.sin(degYZ_Rad) * (pZ - oZ)
        v[2] = oZ + math.sin(degYZ_Rad) * (pY - oY) + math.cos(degYZ_Rad) * (pZ - oZ)
