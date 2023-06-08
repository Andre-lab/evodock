#!/usr/bin/env python
# coding: utf-8

import glob
import random

from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.core.pose import (
    append_pose_to_pose,
    chain_end_res,
)
from pyrosetta.rosetta.core.import_pose import poses_from_files, pose_from_file
from pyrosetta.rosetta.utility import vector1_std_string
from pyrosetta import Pose
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pathlib import Path
from pyrosetta.rosetta.core.pose.symmetry import extract_asymmetric_unit
from scipy.spatial.transform import Rotation as R
from src.utils import get_rotation_euler, get_translation
from src.position_utils import to_rosetta


class FlexbbSwapOperator:
    def __init__(self, config, scfxn, local_search_strategy):
        self.bb_strategy = "library"
        self.scfxn = scfxn
        self.native_fold_tree = scfxn.dock_pose.fold_tree()
        self.config = config
        if self.config.syminfo:
            self.symmetrymover = SetupForSymmetryMover(self.config.syminfo.input_symdef)
        else:
            self.symmetrymover = None
        self.local_search_strategy = local_search_strategy
        if config.syminfo:
            self.list_subunits = []
            if self.config.low_memory_mode:
                self.list_subunits = config.subunit_paths
                # check that we can read all the poses without problems
                for path in self.list_subunits:
                    try:
                        pose_from_file(path)
                    except Exception as e:
                        print(e)
                        exit()
            else:
                filenames_subunits = vector1_std_string()
                for f in config.subunit_paths:
                    filenames_subunits.append(f)
                for pose in poses_from_files(filenames_subunits):
                    self.list_subunits.append(pose) # symmetrize later
                # fixme: check if this is needed. Daniel says that he has problem with the datacache not changning when calculating ZD,
                #  if these are already in the pdb when read.
                self.list_subunits = [Pose(p) for p in self.list_subunits]
            # the relaxed backbones are newer used?
            # self.relaxed_backbones_subunits = [p.clone() for p in self.list_subunits]

        else:
            lst_ligand = glob.glob(config.path_to_ligand_folder)
            lst_receptor = glob.glob(config.path_to_receptor_folder)
            filenames_ligand = vector1_std_string()
            for f in lst_ligand:
                filenames_ligand.append(f)
            filenames_receptor = vector1_std_string()
            for f in lst_receptor:
                filenames_receptor.append(f)
            self.list_ligand = poses_from_files(filenames_ligand)
            self.list_receptor = poses_from_files(filenames_receptor)
            # fixme: check if this is needed. Daniel says that he has problem with the datacache not changning when calculating ZD
            #  if these are already in the pdb when read.
            self.list_ligand = [Pose(p) for p in self.list_ligand]
            self.list_receptor = [Pose(p) for p in self.list_receptor]
            # print("list ligand {}".format(len(self.list_ligand)))
            # print("list receptor {}".format(len(self.list_receptor)))
            # self.relaxed_backbones_ligand = [Pose(p) for p in self.list_ligand]
            # self.relaxed_backbones_receptor = [Pose(p) for p in self.list_receptor]

        # fixme: since this is setting 1 dof make it faster by either working on the translation or the rotation

    def set_symmetric_jump_dof(self, pose, jump_id, dof: int, value):
        """Set the jump with dof to a value"""
        flexible_jump = pose.jump(jump_id)
        rot = get_rotation_euler(flexible_jump)
        trans = get_translation(flexible_jump)
        if dof < 4:
            trans[dof - 1] = value
        else:
            rot[dof - 4] = value
        # convert to Rosetta and set in jump
        rot = R.from_euler("xyz", rot, degrees=True).as_matrix()
        rot, trans = to_rosetta(rot, trans)
        flexible_jump.set_translation(trans)
        flexible_jump.set_rotation(rot)
        pose.set_jump(jump_id, flexible_jump)
        return pose

    def extract_asymmetric_pose(self, pose, syminfo):
        """Extract the asymmetric pose from the symmetric pose and reset the genotype dofs to 0."""
        pose = pose.clone()  # will shadow the pose and pose will not change
        # set all degrees of freedom to 0
        for jump, dofs in zip(syminfo.dof_spec.jump_int, syminfo.dof_spec.dof_int):
            for dof in dofs:
                self.set_symmetric_jump_dof(pose, jump, dof, 0)
        # create a new pose object and fill it with the asymmetric pose
        apose = Pose()
        extract_asymmetric_unit(pose, apose, False)
        return apose

    # def save_bb(self, pose, idx_l, idx_r, idx_s):
    #     """Saves the backbones to file."""
    #     if is_symmetric(pose):
    #         path = self.list_subunits[idx_s]
    #         self.extract_asymmetric_pose(pose, self.config.syminfo).dump_pdb(path)
    #     else:
    #         raise NotImplementedError

    def define_ensemble(self, ind, reference_pose, randomize_bb=True):
        improve_relax = False
        idx_receptor, idx_ligand, idx_subunit = ind.idx_receptor, ind.idx_ligand, ind.idx_subunit
        if randomize_bb is True:
            if is_symmetric(reference_pose):
                idx_subunit = random.randint(0, len(self.list_subunits) - 1)
            else:
                if random.uniform(0, 1) < 0.5:
                    idx_receptor = ind.idx_receptor
                    idx_ligand = random.randint(0, len(self.list_ligand) - 1)
                else:
                    idx_receptor = random.randint(0, len(self.list_receptor) - 1)
                    idx_ligand = ind.idx_ligand

        if is_symmetric(reference_pose):
            join_pose = self.list_subunits[idx_subunit]
            if self.config.low_memory_mode:
                join_pose = pose_from_file(join_pose)
            else:
                # cloning is to make sure we don't symmetrize the stored pose in self.list_subunits
                join_pose = join_pose.clone()
            if not is_symmetric(join_pose):
                self.symmetrymover.apply(join_pose)
            for jump in self.config.syminfo.dof_spec.jump_int:
                join_pose.set_jump(jump, reference_pose.jump(jump))
        else:
            pose_chainA = self.list_receptor[idx_receptor]
            pose_chainB = self.list_ligand[idx_ligand]
            join_pose = self.make_pose_with_chains(reference_pose, pose_chainA, pose_chainB)
        return join_pose, idx_receptor, idx_ligand, idx_subunit, improve_relax

    def make_pose_with_chains(self, reference_pose, pose_chainA, pose_chainB):
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
        join_pose.pdb_info().rebuild_pdb2pose()
        join_pose.fold_tree(self.native_fold_tree)
        join_pose.conformation().detect_disulfides()
        return join_pose

    def apply_bb_strategy(self, ind, pose):
        join_pose, idx_receptor, idx_ligand, idx_subunit, improve_relax = self.define_ensemble(
            ind, pose
        )
        return join_pose, idx_receptor, idx_ligand, idx_subunit, improve_relax
