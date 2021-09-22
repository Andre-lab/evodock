#!/usr/bin/env python
# coding: utf-8


import numpy as np
from pyrosetta import Pose, pose_from_file
from scipy.spatial.transform import Rotation as R

from src.pdb_structure import pdbstructure_from_file

IP_ADDRESS = "10.8.0.6"


def get_pose_from_file(pose_input):
    native_pose = Pose()
    pose_from_file(native_pose, pose_input)
    return native_pose


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
    flexible_jump = dock_pose.jump(dock_pose.num_jump())
    euler_vec = get_rotation_euler(flexible_jump)
    return list(euler_vec) + list(flexible_jump.get_translation())


def set_new_max_translations(scfxn, popul):
    # 1) get all positions
    # 2) set new max_translations
    # 3) convert genotypes to new translation
    all_x = []
    all_y = []
    all_z = []
    all_positions = []
    for ind in popul:
        pose = scfxn.apply_genotype_to_pose(ind.genotype)
        positions = get_position_info(pose)
        all_x.append(positions[3])
        all_y.append(positions[4])
        all_z.append(positions[5])
        all_positions.append(positions)

    max_x = max([val for val in all_x])
    max_x += max_x * 0.1
    max_y = max([val for val in all_y])
    max_y += max_y * 0.1
    max_z = max([val for val in all_z])
    max_z += max_z * 0.1

    min_x = min([val for val in all_x])
    min_x += min_x * 0.1
    min_y = min([val for val in all_y])
    min_y += min_y * 0.1
    min_z = min([val for val in all_z])
    min_z += min_z * 0.1

    # print("prev max_trans ")
    # print(scfxn.converter.max_trans)
    scfxn.converter.max_trans = [max_x, max_y, max_z]
    scfxn.converter.min_trans = [min_x, min_y, min_z]
    scfxn.converter.bounds = scfxn.converter.define_bounds()
    # print("updated max_trans ")
    # print(scfxn.converter.max_trans)
    # print(scfxn.converter.min_trans)

    for i, ind in enumerate(popul):
        # print("=== before ===")
        # print(all_positions[i])
        ind.genotype = scfxn.convert_positions_to_genotype(all_positions[i])
        # print("==== after ===")
        # print(scfxn.convert_genotype(ind.genotype))

    return popul
