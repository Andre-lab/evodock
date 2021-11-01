#!/usr/bin/env python
# coding: utf-8

import glob
import random
from random import sample

from pyrosetta import Pose
from pyrosetta.rosetta.core.import_pose import poses_from_files
from pyrosetta.rosetta.core.pose import (append_pose_to_pose, chain_end_res,
                                         remove_virtual_residues)
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.protocols.docking import (DockMCMProtocol,
                                                 FaDockingSlideIntoContact)
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.utility import vector1_std_string


class LocalSearchStrategy:
    def __init__(self, config, scfxn, dock_pose):
        self.config = config
        self.scfxn = scfxn
        self.dock_pose = dock_pose
        self.packer_option = config.local_search_option
        self.native_fold_tree = scfxn.dock_pose.fold_tree()

        if self.packer_option == "mcm_rosetta":
            mcm_docking = DockMCMProtocol()
            mcm_docking.set_native_pose(scfxn.dock_pose)
            mcm_docking.set_scorefxn(scfxn.scfxn_rosetta)
            mcm_docking.set_rt_min(False)
            mcm_docking.set_sc_min(False)
            mock_pose = Pose()
            mock_pose.assign(scfxn.dock_pose)
            mcm_docking.apply(mock_pose)
            self.docking = mcm_docking
            self.docking.set_task_factory(mcm_docking.task_factory())
            self.docking.set_ignore_default_task(True)

        self.slide_into_contact = FaDockingSlideIntoContact(dock_pose.num_jump())
        if config.docking_type_option == "Unbound":
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
            print("list ligand {}".format(len(self.list_ligand)))
            print("list receptor {}".format(len(self.list_receptor)))
            self.relax = FastRelax(0)
            self.relax.set_scorefxn(self.scfxn.scfxn_rosetta)
            self.relax.max_iter(1)
            self.relaxed_backbones_ligand = [Pose(p) for p in self.list_ligand]
            self.relaxed_backbones_receptor = [Pose(p) for p in self.list_receptor]

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def apply_bb_strategy(self, ind, pose):
        bb_strategy = self.config.bb_strategy
        if bb_strategy == "library":
            join_pose, idx_receptor, idx_ligand = self.define_ensemble(ind, pose)
            return join_pose, idx_receptor, idx_ligand
        if bb_strategy == "relax":
            if random.uniform(0, 1) < 0.5:
                join_pose, idx_receptor, idx_ligand = self.define_relaxedbackbone(
                    ind, pose
                )
            else:
                join_pose, idx_receptor, idx_ligand = self.define_ensemble(ind, pose)
            return join_pose, idx_receptor, idx_ligand
        if bb_strategy == "fixed":
            join_pose = pose
            idx_receptor, idx_ligand = 1, 1
        return join_pose, idx_receptor, idx_ligand

    def add_relaxed_backbones_to_list(self, ind, current_score, reference_pose):
        if current_score < ind.score:
            pose_receptor = Pose(reference_pose, 1, chain_end_res(reference_pose, 1))
            pose_ligand = Pose(
                reference_pose,
                chain_end_res(reference_pose, 1) + 1,
                reference_pose.total_residue(),
            )
            idx_receptor, idx_ligand = ind.idx_receptor, ind.idx_ligand
            idx_receptor = min(idx_receptor, len(self.relaxed_backbones_receptor) - 1)
            idx_ligand = min(idx_ligand, len(self.relaxed_backbones_ligand) - 1)
            self.relaxed_backbones_receptor[idx_receptor] = pose_receptor
            self.relaxed_backbones_ligand[idx_ligand] = pose_ligand

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

    def define_ensemble(self, ind, reference_pose, randomize_bb=True):
        if randomize_bb is True:
            if random.uniform(0, 1) < 0.5:
                idx_receptor = ind.idx_receptor
                idx_ligand = random.randint(0, len(self.list_ligand) - 1)
            else:
                idx_receptor = random.randint(0, len(self.list_receptor) - 1)
                idx_ligand = ind.idx_ligand
        else:
            idx_receptor, idx_ligand = ind.idx_receptor, ind.idx_ligand

        idx_receptor = min(idx_receptor, len(self.list_receptor) - 1)
        idx_ligand = min(idx_ligand, len(self.list_ligand) - 1)
        pose_chainA = self.list_receptor[idx_receptor]
        pose_chainB = self.list_ligand[idx_ligand]
        join_pose = self.make_pose_with_chains(reference_pose, pose_chainA, pose_chainB)

        return join_pose, idx_receptor, idx_ligand

    def define_relaxedbackbone(self, ind, pose):
        idx_receptor, idx_ligand = ind.idx_receptor, ind.idx_ligand
        idx_receptor = min(idx_receptor, len(self.relaxed_backbones_receptor) - 1)
        idx_ligand = min(idx_ligand, len(self.relaxed_backbones_ligand) - 1)

        pose_chainA = self.relaxed_backbones_receptor[idx_receptor]
        pose_chainB = self.relaxed_backbones_ligand[idx_ligand]
        join_pose = self.make_pose_with_chains(pose, pose_chainA, pose_chainB)
        return join_pose, idx_receptor, idx_ligand

    def apply_bound_docking(self, ind, local_search=True):
        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
        before = self.energy_score(pose)
        if local_search and self.packer_option != "None":
            self.slide_into_contact.apply(pose)
            self.docking.apply(pose)
            after = self.energy_score(pose)
        else:
            after = before

        return_data = {
            "pose": pose,
            "before": before,
            "after": after,
            "idx_ligand": 1,
            "idx_receptor": 1,
        }
        return return_data

    def apply_unbound_docking(self, ind, local_search=True):
        relax_prob = self.config.relax_prob
        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
        join_pose, idx_receptor, idx_ligand = self.apply_bb_strategy(ind, pose)
        before = self.energy_score(join_pose)
        if local_search and self.packer_option != "None":
            apply_relax = relax_prob > random.uniform(0, 1)
            if not apply_relax:
                self.slide_into_contact.apply(join_pose)
                self.docking.apply(join_pose)
                after = self.energy_score(join_pose)
            else:
                self.slide_into_contact.apply(join_pose)
                self.docking.apply(join_pose)
                self.relax.apply(join_pose)
                remove_virtual_residues(join_pose)
                join_pose.fold_tree(self.native_fold_tree)
                after = self.energy_score(join_pose)
                self.add_relaxed_backbones_to_list(ind, after, join_pose)

        else:
            after = before
        return_data = {
            "pose": join_pose,
            "before": before,
            "after": after,
            "idx_ligand": idx_ligand,
            "idx_receptor": idx_receptor,
        }
        return return_data

    def apply(self, ind, local_search=True):
        if self.config.docking_type_option == "Unbound":
            return self.apply_unbound_docking(ind, local_search)
        else:
            return self.apply_bound_docking(ind, local_search)
