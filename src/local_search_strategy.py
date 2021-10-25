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
from pyrosetta.rosetta.protocols.docking import (DockingSlideIntoContact,
                                                 DockMCMProtocol)
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

        self.slide_into_contact = DockingSlideIntoContact(dock_pose.num_jump())
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
            print("list ligand {}".format(len(self.list_ligand)))
            print("list receptor {}".format(len(self.list_receptor)))
            self.relax = FastRelax(0)
            self.relax.set_scorefxn(self.scfxn.scfxn_rosetta)
            self.relax.max_iter(1)
            self.relaxed_backbones = [self.scfxn.dock_pose]

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def apply_bb_strategy(self, ind, pose):
        bb_strategy = self.config.bb_strategy
        if bb_strategy == "library":
            join_pose, idx_receptor, idx_ligand = self.define_ensemble(ind, pose)
            return join_pose, idx_receptor, idx_ligand
        if bb_strategy == "relax":
            if self.config.relax_prob > random.uniform(0, 1):
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

    def define_ensemble(self, ind, reference_pose, randomize_bb=True):
        if randomize_bb is True:
            if ind.idx_receptor == ind.idx_ligand:
                idx_ligand = random.randint(1, len(self.list_ligand))
                idx_receptor = random.randint(1, len(self.list_receptor))
            else:
                if random.uniform(0, 1) < 0.5:
                    idx_receptor = ind.idx_receptor
                    idx_ligand = random.randint(1, len(self.list_ligand))
                else:
                    idx_receptor = random.randint(1, len(self.list_receptor))
                    idx_ligand = ind.idx_ligand
        else:
            idx_receptor, idx_ligand = ind.idx_receptor, ind.idx_ligand

        pose_chainA = self.list_receptor[idx_receptor]
        pose_chainB = self.list_ligand[idx_ligand]

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

        return join_pose, idx_receptor, idx_ligand

    def define_relaxedbackbone(self, ind, pose):
        print("relax backbone")
        idx = random.randint(0, len(self.relaxed_backbones) - 1)
        if idx > len(self.relaxed_backbones):
            idx = random.randint(0, len(self.relaxed_backbones) - 1)

        join_pose = self.relaxed_backbones[idx]
        idx_receptor, idx_ligand = idx, idx

        calpha_superimpose_pose(join_pose, pose)
        join_pose.fold_tree(self.native_fold_tree)
        join_pose.conformation().detect_disulfides()
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
            else:
                self.relax.apply(join_pose)
                self.relaxed_backbones.append(join_pose)
                print("backbone")
                if len(self.relaxed_backbones) > 100:
                    self.relaxed_backbones = sample(self.relaxed_backbones, 100)
                idx_ligand, idx_receptor = (
                    len(self.relaxed_backbones),
                    len(self.relaxed_backbones),
                )
                remove_virtual_residues(join_pose)
                join_pose.fold_tree(self.native_fold_tree)
            after = self.energy_score(join_pose)
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
