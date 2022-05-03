#!/usr/bin/env python
# coding: utf-8


import numpy as np
from pyrosetta import Pose, pose_from_file, Vector1
from scipy.spatial.transform import Rotation as R
from src.individual import Individual
from src.pdb_structure import pdbstructure_from_file
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_jump_num, jump_num_sym_dof
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pose import chain_end_res
from pyrosetta.rosetta.protocols.docking import setup_foldtree

IP_ADDRESS = "10.8.0.6"

# symmetric maps from dof name to dof int
strtodofint = {"x": 1, "y": 2, "z": 3, "angle_x": 4, "angle_y": 5, "angle_z": 6}
inttodofstr =  {v: k for k, v in strtodofint.items()}

def map_jump_and_dofs_to_int(pose: Pose, syminfo: dict) -> None:
    """Map jump strings to jump ints and dof str to dof ints."""
    syminfo["jumps_int"] = [sym_dof_jump_num(pose, jump_str) for jump_str in syminfo.get("jumps_str")]
    syminfo["dofs_int"] = [[strtodofint.get(dof_str) for dof_str in dofs_str] for dofs_str in syminfo.get("dofs_str")]

def get_symmetric_genotype_str(pose: Pose) -> str:
    """Gets which symdofs are present in the genotype and in the order as they are in the genotype."""
    symdofs = pose.conformation().Symmetry_Info().get_dofs()
    symmetric_genotype = []
    for jump_id, symdof in symdofs.items():
        jump_str = jump_num_sym_dof(pose, jump_id)
        for dof in range(1, 7):
            if symdof.allow_dof(dof):
                symmetric_genotype.append(f"{jump_str}:{inttodofstr[dof]}")
    return " ".join(symmetric_genotype)

def get_starting_poses(pose_input, native_input, config):
    native = Pose()
    pose_from_file(native, native_input)
    pose = Pose()
    pose_from_file(pose, pose_input)
    native.conformation().detect_disulfides()
    pose.conformation().detect_disulfides()
    if config.syminfo:
        SetupForSymmetryMover(config.syminfo.get("input_symdef")).apply(pose)
        SetupForSymmetryMover(config.syminfo.get("native_symdef")).apply(native)
        map_jump_and_dofs_to_int(pose, config.syminfo)
    else:
        mres = chain_end_res(pose, 1)
        ft = FoldTree()
        ft.add_edge(1, mres, -1)
        ft.add_edge(1, mres + 1, 1)
        ft.add_edge(mres + 1, pose.total_residue(), -1)
        pose.fold_tree(ft)
    if config.pmm:
        pose.pdb_info().name("input_pose")
        native.pdb_info().name("native_pose")
        config.pmm.apply(pose)
        config.pmm.apply(native)
    return pose, native

def get_pose_from_file(pose_input):
    pose = Pose()
    pose_from_file(pose, pose_input)
    pose.conformation().detect_disulfides()
    return pose
# >>>>>>> main

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

def get_translation(flexible_jump):
    return np.asarray(flexible_jump.rt().get_translation())

def get_position_info(dock_pose):
    if is_symmetric(dock_pose):
        native_val = []
        dofs = dock_pose.conformation().Symmetry_Info().get_dofs()
        for jump_id, symdof in dofs.items():
            flexible_jump = dock_pose.jump(jump_id)
            rot = get_rotation_euler(flexible_jump)
            trans = get_translation(flexible_jump)
            for dof in range(1, 7):
                if symdof.allow_dof(dof):
                    if dof < 4:
                        native_val.append(trans[dof - 1])
                    else:
                        native_val.append(rot[dof - 4])
    else:
        flexible_jump = dock_pose.jump(1)
        euler_vec = get_rotation_euler(flexible_jump)
        native_val = list(euler_vec) + list(flexible_jump.get_translation())
    return native_val


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


def make_trial(idx, genotype, ligand=1, receptor=1, subunit=1):
    ind = Individual(idx, genotype)
    ind.idx = idx
    ind.genotype = genotype
    ind.score = 1000
    ind.idx_subunit = subunit
    ind.idx_ligand = ligand
    ind.idx_receptor = receptor
    ind.rmsd = 0
    ind.i_sc = 0
    ind.irms = 0
    return ind
