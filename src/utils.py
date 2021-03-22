#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.spatial.transform import Rotation as R

from src.pdb_structure import pdbstructure_from_file

IP_ADDRESS = "10.8.0.6"


# compute an axis-aligned bounding box for the given pdb structure
def xyz_limits_for_pdb(pdb):
    first = True
    count = 0
    for chain in pdb.chains:
        for res in chain.residues:
            for atom in res.atoms:
                count += 1
                if first:
                    # print "first", count
                    first = False
                    lower_xyz = np.array(atom.xyz).tolist()
                    upper_xyz = np.array(atom.xyz).tolist()
                else:
                    lower_xyz = np.minimum(lower_xyz, np.array(atom.xyz).tolist())
                    upper_xyz = np.maximum(upper_xyz, np.array(atom.xyz).tolist())
    # print "xyz from", count, "atoms"
    return lower_xyz, upper_xyz


def get_translation_max(input_pdb):
    pdb = pdbstructure_from_file(input_pdb)
    lower_xyz, upper_xyz = xyz_limits_for_pdb(pdb)
    max_translation = max(np.maximum(abs(lower_xyz), abs(upper_xyz)))
    return max_translation + 10


def convert_range(OldValue, old_range, new_range):
    OldMin = old_range[0]
    OldMax = old_range[1]
    NewMin = new_range[0]
    NewMax = new_range[1]
    NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
    return NewValue


def convert_translation(OldValue, trans_max_magnitude):
    OldMin = -1
    OldMax = 1
    NewMin = -1 * trans_max_magnitude
    NewMax = trans_max_magnitude
    NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
    return NewValue


def convert_gamma(OldValue):
    return ((OldValue - (-1)) * (180)) / (1 - (-1))


def random_individual(max_value=180):
    ind = []
    for i in range(6):
        if i < 1:
            # ind.append(0)
            ind.append(np.random.random_sample() * max_value)
        else:
            # ind.append(0)
            ind.append(np.random.random_sample() * 1)
    return ind


def get_rotation_euler(flexible_jump):
    rot_matrix = R.from_matrix(np.asarray(flexible_jump.get_rotation()))
    vec = rot_matrix.as_euler("xyz", degrees=True)
    return vec


def get_position_info(dock_pose):
    flexible_jump = dock_pose.jump(1)
    euler_vec = get_rotation_euler(flexible_jump)
    return list(euler_vec) + list(flexible_jump.get_translation())
