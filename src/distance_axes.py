#! /usr/bin/env python

"""
Computes principal axes from a PDB file.

Produces also a .pml script for a nice rendering with PyMOL.
"""

import numpy as np
from pyrosetta.rosetta.core.pose import get_resnums_for_chain_id

__author__ = "Pierre Poulain"
__credits__ = ["Justine Guegan", "Edithe Selwa", "Steven C. Howell"]
__license__ = "GPL"


# scale factor to enhance the length of axis in Pymol
scale_factor = 20


axes = {"x": 0, "y": 1, "z": 2}


def read_pdb_xyz(pose, chain, ax):
    distances = []
    residues = list(get_resnums_for_chain_id(pose, chain))
    for i in range(0, len(residues) - 1):
        for j in range(i, len(residues)):
            x1 = pose.residue(residues[i]).xyz("CA")[axes[ax]]
            x2 = pose.residue(residues[j]).xyz("CA")[axes[ax]]
            distances.append(abs(x1 - x2))
    return max(distances)


def calculate_local_coordinates(pose):
    jump_num = 1
    max_trans_amount = 3
    flexible_jump = pose.jump(jump_num)
    translation = np.asarray(flexible_jump.get_translation())
    max_trans = [translation[i] + max_trans_amount for i in range(len(translation))]
    min_trans = [translation[i] - max_trans_amount for i in range(len(translation))]
    return max_trans, min_trans


def calculate_max_coordiantes(pose):
    max_x_1 = read_pdb_xyz(pose, 1, "x")
    max_y_1 = read_pdb_xyz(pose, 1, "y")
    max_z_1 = read_pdb_xyz(pose, 1, "z")

    max_x_2 = read_pdb_xyz(pose, 2, "x")
    max_y_2 = read_pdb_xyz(pose, 2, "y")
    max_z_2 = read_pdb_xyz(pose, 2, "z")

    min_x = max_x_1 / 2
    min_y = max_y_1 / 2
    min_z = max_z_1 / 2
    # max_x = (max_x_1 / 2) + (max_x_2 / 2) + 10
    # max_y = (max_y_1 / 2) + (max_y_2 / 2) + 10
    # max_z = (max_z_1 / 2) + (max_z_2 / 2) + 10

    max_x = (max_x_1 / 2) + (max_x_2)
    max_y = (max_y_1 / 2) + (max_y_2)
    max_z = (max_z_1 / 2) + (max_z_2)

    # print([min_x, min_y, min_z])
    max_all = max([max_x, max_y, max_z])
    min_all = max_all * -1
    return [max_all, max_all, max_all], [min_all, min_all, min_all]
