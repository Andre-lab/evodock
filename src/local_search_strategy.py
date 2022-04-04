#!/usr/bin/env python
# coding: utf-8

import random

from src.flexbb_swap_operator import FlexbbSwapOperator
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from pyrosetta import Pose
from pyrosetta.rosetta.protocols.docking import (
    DockMCMProtocol,
    FaDockingSlideIntoContact,
)


from pyrosetta import Vector1
from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent,
    InitializeFromCommandline,
    NoRepackDisulfides,
    RestrictToRepacking,
)
from pyrosetta.rosetta.protocols.simple_task_operations import RestrictToInterface


class LocalSearchStrategy:
    def __init__(self, config, scfxn, dock_pose):
        self.config = config
        # self.pymover = PyMOLMover(address="10.8.0.22", port=65000, max_packet_size=1400)
        self.scfxn = scfxn
        self.dock_pose = dock_pose
        self.packer_option = config.local_search_option
        self.native_fold_tree = scfxn.dock_pose.fold_tree()
        if self.packer_option == "sidechains":
            local_tf = TaskFactory()
            local_tf.push_back(InitializeFromCommandline())
            local_tf.push_back(IncludeCurrent())
            local_tf.push_back(RestrictToRepacking())
            local_tf.push_back(NoRepackDisulfides())
            extrarot = operation.ExtraRotamersGeneric()
            extrarot.ex1(True)
            extrarot.ex2aro(True)
            local_tf.push_back(extrarot)
            restrict = RestrictToInterface()
            restrict.set_movable_jumps(Vector1([1]))
            local_tf.push_back(restrict)
            mcm_docking = DockMCMProtocol()
            mcm_docking.set_native_pose(scfxn.dock_pose)
            mcm_docking.set_scorefxn(scfxn.scfxn_rosetta)
            mcm_docking.set_rt_min(False)
            mcm_docking.set_sc_min(False)
            mock_pose = Pose()
            mock_pose.assign(scfxn.dock_pose)
            mcm_docking.apply(mock_pose)
            self.docking = mcm_docking
            self.docking.set_task_factory(local_tf)
            self.docking.set_ignore_default_task(True)
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
        if self.config.docking_type_option == "Unbound":
            self.swap_operator = FlexbbSwapOperator(config, scfxn, None)

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def apply_bb_strategy(self, ind, pose):
        bb_strategy = self.config.bb_strategy
        if bb_strategy == "library":
            return self.swap_operator.apply_bb_strategy(ind, pose)
        if bb_strategy == "popul_library":
            if random.uniform(0, 1) < 0.3:
                idx_receptor, idx_ligand = ind.idx_receptor, ind.idx_ligand
                pose_chainA = self.swap_operator.list_receptor[idx_receptor]
                pose_chainB = self.swap_operator.list_ligand[idx_ligand]
                join_pose = self.swap_operator.make_pose_with_chains(
                    pose, pose_chainA, pose_chainB
                )
                return join_pose, idx_receptor, idx_ligand

        join_pose = pose
        idx_receptor, idx_ligand = ind.idx_receptor, idx_ligand
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
            "idx_ligand": ind.idx_ligand,
            "idx_receptor": ind.idx_receptor,
        }
        return return_data

    def apply_unbound_docking(self, ind, local_search=True):
        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
        # pose.pdb_info().name("apply_gen_" + str(ind.idx))
        # self.pymover.apply(pose)
        join_pose, idx_receptor, idx_ligand = self.apply_bb_strategy(ind, pose)
        # join_pose.pdb_info().name("apply_bb_" + str(ind.idx))
        # print("apply_bb_" + str(ind.idx) + " : " + str(self.energy_score(join_pose)))
        # self.pymover.apply(join_pose)
        before = self.energy_score(join_pose)
        if local_search and self.packer_option != "None":
            self.slide_into_contact.apply(join_pose)
            self.docking.apply(join_pose)
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
