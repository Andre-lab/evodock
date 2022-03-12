#!usr/bin/env python

from __future__ import print_function

import glob
import time

import numpy as np
import pyrosetta.rosetta as rosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.pose import addVirtualResAsRoot, append_pose_to_pose
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from pyrosetta.rosetta.protocols.rigid import RigidBodySpinMover
from pyrosetta.rosetta.protocols.toolbox.rigid_body import create_euler_rotation
from scipy.spatial import distance


def print_rotation_translation(R, t):
    with open("align_info.txt", "w") as file_object:
        for j in R:
            for i in j:
                file_object.write(str(i).strip() + "\n")
        for i in t:
            file_object.write(str(i) + "\n")


def to_rosetta(rotation, translation):
    mros = rosetta.numeric.xyzMatrix_double_t(0)
    vec = rosetta.numeric.xyzVector_double_t(0)

    for i in [0, 1, 2]:
        vec.assign(rotation[i][0], rotation[i][1], rotation[i][2])
        mros.row(i + 1, vec)

    trans_vector = rosetta.numeric.xyzVector_double_t(0)
    trans_vector.assign(translation[0], translation[1], translation[2])
    return mros, trans_vector


def read_rotation_translation():
    with open("align_info.txt", "r") as file_object:
        lines = file_object.readlines()
        rotation = np.array(lines[:9], dtype="float64").reshape(3, 3)
        translation = np.array(lines[9:], dtype="float64").reshape(1, 3)[0]

    return to_rosetta(rotation, translation)


def build_rotation_and_translation():
    mros = rosetta.numeric.xyzMatrix_double_t(0)
    vec = rosetta.numeric.xyzVector_double_t(0)
    vec.assign(
        0.9523832594972914,
        -0.1968083849638013,
        -0.2328789098163566,
    )
    mros.row(0, vec)
    vec.assign(-0.03736108525272122, -0.8333504234478073, 0.5514809344375380)
    mros.row(1, vec)
    vec.assign(-0.3026058101525034, -0.5165206010870249, -0.8010219674357533)
    mros.row(2, vec)

    trans_vector = rosetta.numeric.xyzVector_double_t(0)
    trans_vector.assign(-11.28210422438503, -11.92831975239535, 9.604356732570816)
    # trans_vector.assign(0,0,0)

    return mros, trans_vector


def get_translation(x, y, z):
    axis1 = rosetta.numeric.xyzVector_double_t(0)
    axis1.assign(x, y, z)
    return axis1


def build_axis():
    axis1 = rosetta.numeric.xyzVector_double_t(0)
    axis1.assign(1, 0, 0)
    axis2 = rosetta.numeric.xyzVector_double_t(0)
    axis2.assign(0, 1, 0)
    axis3 = rosetta.numeric.xyzVector_double_t(0)
    axis3.assign(0, 0, 1)
    return axis1, axis2, axis3


def get_anchor_coordinates_from_pose(residue):
    bb_atoms = ["N", "CA", "C"]
    coords = []
    for atom in bb_atoms:
        coords.append([residue.xyz(atom).x, residue.xyz(atom).y, residue.xyz(atom).z])
    return np.mat(coords)


def build_start_mros_trans():
    mros = rosetta.numeric.xyzMatrix_double_t(0)
    vec = rosetta.numeric.xyzVector_double_t(0)
    vec.assign(1, 0, 0)
    mros.row(0, vec)
    vec.assign(0, 1, 0)
    mros.row(1, vec)
    vec.assign(0, 0, 1)
    mros.row(2, vec)

    trans_vector = rosetta.numeric.xyzVector_double_t(0)
    trans_vector.assign(0, 0, 0)
    # trans_vector.assign(0,0,0)

    return mros, trans_vector


def print_rotation(pose1, pose2, i):
    join_pose = Pose()
    join_pose.assign(pose1)
    append_pose_to_pose(join_pose, pose2, True)
    join_pose.dump_pdb("Sample_" + str(i) + ".pdb")


def show_pose(pose1, pose2):
    # pymover = PyMOLMover()
    join_pose = Pose()
    join_pose.assign(pose1)
    append_pose_to_pose(join_pose, pose2, True)

    jump_num = 1
    print("rotation")
    print(join_pose.jump(jump_num).get_rotation())  # rotation matrix
    print("translation")
    print(join_pose.jump(jump_num).get_translation())  # translation vector
    # pymover.apply(join_pose)


def get_ax1(i):
    axis1 = rosetta.numeric.xyzVector_double_t(0)
    axis1.assign(i, 0, 0)
    return axis1
