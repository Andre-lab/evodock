#!/usr/bin/env python
# coding: utf-8
import numpy as np
from pyrosetta import Pose, pose_from_file, Vector1
from pyrosetta.rosetta.basic.datacache import CacheableStringMap
from pyrosetta.rosetta.core.pose.datacache import CacheableDataType
from scipy.spatial.transform import Rotation as R
from src.individual import Individual
from src.pdb_structure import pdbstructure_from_file
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.pose import chain_end_res
from pathlib import Path
from symmetryhandler.reference_kinematics import set_jumpdof_str_str
from cubicsym.utilities import get_jumpidentifier
from cubicsym.utilities import get_base_from_pose, add_id_to_pose_w_base
from cubicsym.actors.symdefswapper import SymDefSwapper
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose

def set_init_x(pose, config):
    """Sets the initial x for the input pose if an x_transfile is given in the config file"""
    xf = config.syminfo.x_transfile
    if xf is not None:
        init_x = xf[xf["model"] == Path(config.pose_input).name]["x_trans"].values[0]
        set_jumpdof_str_str(pose, f"JUMP{get_jumpidentifier(pose)}fold111", "x", init_x)

def initialize_starting_poses(config):
    """Initialize the starting poses."""
    native = None
    if config.native_input is not None:
        native = Pose()
        pose_from_file(native, config.native_input)
        native.conformation().detect_disulfides()
    # if flexbb we have to merge the receptor and ligand together
    if isinstance(config.pose_input, list):
        receptor = pose_from_file(config.pose_input[0])
        ligand = pose_from_file(config.pose_input[1])
        if config.template is not None:
            input_pose = align_two_chains_to_pose(pose_from_file(config.template), receptor, ligand)
        else:
            append_pose_to_pose(receptor, ligand, True)
            input_pose = receptor.clone()
    else:
        input_pose = Pose()
        pose_from_file(input_pose, config.pose_input)
    input_pose.conformation().detect_disulfides()
    native_symmetric = None
    if config.syminfo:
        # native suff
        SetupForSymmetryMover(config.syminfo.input_symdef).apply(input_pose)
        # symmetrize the native so that it is the same base
        if config.syminfo.native_symmetric_input is not None:
            native_symmetric = pose_from_file(config.syminfo.native_symmetric_input)
            native_symmetric.conformation().detect_disulfides()
            SetupForSymmetryMover(config.syminfo.native_symdef).apply(native_symmetric)
            base = get_base_from_pose(input_pose)
            if base != "HF":
                sds = SymDefSwapper(native_symmetric, config.syminfo.native_symdef)
                if base == "3F":
                    native_symmetric = sds.create_3fold_pose(native_symmetric)
                elif base == "2F":
                    native_symmetric = sds.create_2fold_pose(native_symmetric)
                else:
                    raise ValueError
        set_init_x(input_pose, config)
        config.syminfo.store_info_from_pose(input_pose) # setup cubic boundaries
    else:
        mres = chain_end_res(input_pose, 1)
        ft = FoldTree()
        ft.add_edge(1, mres, -1)
        ft.add_edge(1, mres + 1, 1)
        ft.add_edge(mres + 1, input_pose.total_residue(), -1)
        input_pose.fold_tree(ft)
    if config.pmm:
        input_pose.pdb_info().name("input_pose")
        if config.native_input is not None:
            native.pdb_info().name("native_pose")
        config.pmm.apply(input_pose)
        config.pmm.apply(native)
    return input_pose, native, native_symmetric

def get_pose_from_file(pose_input):
    pose = Pose()
    pose_from_file(pose, pose_input)
    pose.conformation().detect_disulfides()
    return pose

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


def get_translation_max(pose):
    # we just do this on the pose instead
    for resi in range(1, pose.size() + 1):
        if resi == 1:
            for atom in pose.residue(resi).atoms():
                val = np.array(atom.xyz()).tolist()
                lower_xyz = val
                upper_xyz = val
        else:
            for atom in pose.residue(resi).atoms():
                val = np.array(atom.xyz()).tolist()
                lower_xyz = np.minimum(lower_xyz, val)
                upper_xyz = np.maximum(upper_xyz, val)
    # pdb = pdbstructure_from_file(input_pdb)
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

def get_position_info(dock_pose, syminfo=None):
    if is_symmetric(dock_pose):
        assert syminfo is not None, "syminfo should be defined if pose is symmetric!"
        # map the symmetrical info in the right order as specified in syminfo
        return syminfo.get_position_info(dock_pose)
    else:
        flexible_jump = dock_pose.jump(1)
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


def make_trial(idx, genotype, ligand=1, receptor=1, subunit=1, receptor_name="", ligand_name="", subunit_name="", flipped=None, fixed=None):
    return Individual(idx, genotype, score=1000, idx_ligand=ligand, idx_receptor=receptor, idx_subunit=subunit, rmsd=0, i_sc=0, irms=0,
        receptor_name=receptor_name, ligand_name=ligand_name, subunit_name=subunit_name, flipped=flipped, fixed=fixed)

def align_two_chains_to_pose(reference_pose, pose_chainA, pose_chainB):
    """Aligs pose_chainA and pose_chainB onto the chain A and B of the reference_pose respectively."""
    pose_receptor = Pose(reference_pose, 1, chain_end_res(reference_pose, 1))
    pose_ligand = Pose(
        reference_pose,
        chain_end_res(reference_pose, 1) + 1,
        reference_pose.total_residue(),
    )
    calpha_superimpose_pose(pose_chainA, pose_receptor)
    calpha_superimpose_pose(pose_chainB, pose_ligand)
    join_pose = Pose()
    join_pose.assign(pose_chainA)
    append_pose_to_pose(join_pose, pose_chainB, True)
    return join_pose