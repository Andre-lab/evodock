#!/usr/bin/env python
# coding: utf-8

import glob
import random


from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.core.pose import (
    append_pose_to_pose,
    chain_end_res,
    remove_virtual_residues,
)
from pyrosetta.rosetta.core.import_pose import poses_from_files
from pyrosetta.rosetta.utility import vector1_std_string
from pyrosetta import Pose


class FlexbbSwapOperator:
    def __init__(self, config, scfxn, local_search):
        self.bb_strategy = "library"
        self.scfxn = scfxn
        self.native_fold_tree = scfxn.dock_pose.fold_tree()

        self.local_search = local_search
        lst_ligand = glob.glob(config.path_ligands)
        lst_receptor = glob.glob(config.path_receptors)
        filenames_ligand = vector1_std_string()
        for f in lst_ligand:
            filenames_ligand.append(f)
        filenames_receptor = vector1_std_string()
        for f in lst_receptor:
            filenames_receptor.append(f)
        self.list_ligand = poses_from_files(filenames_ligand)
        self.list_receptor = poses_from_files(filenames_receptor)
        self.list_ligand = [Pose(p) for p in self.list_ligand]
        self.list_receptor = [Pose(p) for p in self.list_receptor]
        # print("list ligand {}".format(len(self.list_ligand)))
        # print("list receptor {}".format(len(self.list_receptor)))
        self.relaxed_backbones_ligand = [Pose(p) for p in self.list_ligand]
        self.relaxed_backbones_receptor = [Pose(p) for p in self.list_receptor]
        self.relax = FastRelax(1)
        self.relax.set_scorefxn(self.scfxn.scfxn_rosetta)

    def add_relaxed_backbones_to_list(self, ind, current_score, reference_pose):
        pose_receptor = Pose(reference_pose, 1, chain_end_res(reference_pose, 1))
        pose_ligand = Pose(
            reference_pose,
            chain_end_res(reference_pose, 1) + 1,
            reference_pose.total_residue(),
        )
        idx_receptor, idx_ligand = ind.idx_receptor, ind.idx_ligand
        self.list_receptor[idx_receptor] = pose_receptor
        self.list_ligand[idx_ligand] = pose_ligand

    def define_ensemble(self, ind, reference_pose, randomize_bb=True):
        improve_relax = False
        if randomize_bb is True:
            if random.uniform(0, 1) < 0.5:
                idx_receptor = ind.idx_receptor
                idx_ligand = random.randint(0, len(self.list_ligand) - 1)
            else:
                idx_receptor = random.randint(0, len(self.list_receptor) - 1)
                idx_ligand = ind.idx_ligand
        else:
            idx_receptor, idx_ligand = ind.idx_receptor, ind.idx_ligand

        pose_chainA = self.list_receptor[idx_receptor]
        pose_chainB = self.list_ligand[idx_ligand]
        join_pose = self.make_pose_with_chains(reference_pose, pose_chainA, pose_chainB)

        if random.uniform(0, 1) < 0.5:
            self.local_search.local_search_strategy.slide_into_contact.apply(join_pose)
            self.local_search.local_search_strategy.docking.apply(join_pose)
            before = self.scfxn.scfxn_rosetta(join_pose)
            # join_pose.dump_scored_pdb("before.pdb", self.scfxn.scfxn_rosetta)
            self.relax.apply(join_pose)
            remove_virtual_residues(join_pose)
            join_pose.fold_tree(self.native_fold_tree)
            after = self.scfxn.scfxn_rosetta(join_pose)
            # join_pose.dump_scored_pdb("after.pdb", self.scfxn.scfxn_rosetta)
            print("before {} after {}".format(before, after))
            if after < before:
                self.add_relaxed_backbones_to_list(ind, after, join_pose)
                improve_relax = True

        return join_pose, idx_receptor, idx_ligand, improve_relax

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
        join_pose, idx_receptor, idx_ligand, improve_relax = self.define_ensemble(
            ind, pose
        )
        return join_pose, idx_receptor, idx_ligand, improve_relax
