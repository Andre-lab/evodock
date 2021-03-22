#! /usr/bin/env python

"""
Computes principal axes from a PDB file.

Produces also a .pml script for a nice rendering with PyMOL.
"""

import os.path
import sys

import numpy

__author__ = "Pierre Poulain"
__credits__ = ["Justine Guegan", "Edithe Selwa", "Steven C. Howell"]
__license__ = "GPL"


# scale factor to enhance the length of axis in Pymol
scale_factor = 20


def read_pdb_xyz(pose):
    xyz = []
    for i in range(1, pose.total_residue() + 1):
        position = pose.residue(i).xyz("CA")
        xyz.append(position)
    return xyz


def calculate_max_coordiantes(pose):
    # --------------------------------------------------------------------------
    # compute principal axes
    # --------------------------------------------------------------------------
    # read pdb
    xyz = read_pdb_xyz(pose)
    # print("%d CA atomes found in %s" % (len(xyz), pdb_name))

    # create coordinates array
    coord = numpy.array(xyz, float)

    # compute geometric center
    center = numpy.mean(coord, 0)
    # print("Coordinates of the geometric center:\n", center)

    # center with geometric center
    coord = coord - center

    # compute principal axis matrix
    inertia = numpy.dot(coord.transpose(), coord)
    e_values, e_vectors = numpy.linalg.eig(inertia)
    # warning eigen values are not necessary ordered!
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html
    # print("(Unordered) eigen values:")
    # print(e_values)
    # print("(Unordered) eigen vectors:")
    # print(e_vectors)

    # --------------------------------------------------------------------------
    # order eigen values (and eigen vectors)
    #
    # axis1 is the principal axis with the biggest eigen value (eval1)
    # axis2 is the principal axis with the second biggest eigen value (eval2)
    # axis3 is the principal axis with the smallest eigen value (eval3)
    # --------------------------------------------------------------------------
    order = numpy.argsort(e_values)
    eval3, eval2, eval1 = e_values[order]
    axis3, axis2, axis1 = e_vectors[:, order].transpose()
    # print("Inertia axis are now ordered !")

    # --------------------------------------------------------------------------
    # center axes to the geometric center of the molecule
    # and rescale them by order of eigen values
    # --------------------------------------------------------------------------
    # the large vector is the first principal axis
    point1 = 3 * scale_factor * axis1 + center
    # the medium vector is the second principal axis
    point2 = 2 * scale_factor * axis2 + center
    # the small vector is the third principal axis
    point3 = 1 * scale_factor * axis3 + center

    max_coordinate_x = max([point1[0], point2[0], point3[0]])
    max_coordinate_y = max([point1[1], point2[1], point3[1]])
    max_coordinate_z = max([point1[2], point2[2], point3[2]])
    return [max_coordinate_x, max_coordinate_y, max_coordinate_z]
