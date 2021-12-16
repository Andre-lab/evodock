#!/usr/bin/env python
# coding: utf-8

import glob
import random
from random import sample

from src.flexbb_swap_operator import FlexbbSwapOperator
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from pyrosetta import Pose
from pyrosetta.rosetta.core.import_pose import poses_from_files
from pyrosetta.rosetta.core.pose import (append_pose_to_pose, chain_end_res,
                                         remove_virtual_residues)
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.protocols.docking import (DockMCMProtocol,
                                                 FaDockingSlideIntoContact)
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.utility import vector1_std_string


from pyrosetta import Pose, Vector1, standard_packer_task
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import PackerTask, TaskFactory, operation
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent, InitializeFromCommandline, NoRepackDisulfides,
    RestrictToRepacking)
from pyrosetta.rosetta.protocols.docking import DockMCMProtocol
from pyrosetta.rosetta.protocols.minimization_packing import (
    MinMover, PackRotamersMover, RotamerTrialsMover)
from pyrosetta.rosetta.protocols.moves import (JumpOutMover, MonteCarlo,
                                               SequenceMover, TrialMover)
from pyrosetta.rosetta.protocols.rigid import (RigidBodyPerturbMover,
                                               partner_downstream)
from pyrosetta.rosetta.protocols.simple_task_operations import \
    RestrictToInterface



class LocalSearchStrategy:
    def __init__(self, config, scfxn, dock_pose):
        self.config = config
        self.pymover = PyMOLMover(address="10.8.0.22", port=65000, max_packet_size=1400)
        self.scfxn = scfxn
        self.dock_pose = dock_pose
        self.packer_option = config.local_search_option
        self.native_fold_tree = scfxn.dock_pose.fold_tree()
        self.relax_prob = self.config.relax_prob
        if self.config.bb_strategy == "only_relax":
            self.relax_prob = 1.0
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
        if self.config.docking_type_option == "Unbound" and self.config.bb_strategy != "relax_best":
            self.swap_operator = FlexbbSwapOperator(config, scfxn, None)
            self.relax = FastRelax(0)
            self.relax.set_scorefxn(self.scfxn.scfxn_rosetta)
            self.relax.max_iter(1)

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def apply_bb_strategy(self, ind, pose):
        bb_strategy = self.config.bb_strategy
        if bb_strategy == "library":
            return self.swap_operator.apply_bb_strategy(ind, pose)
        if bb_strategy == "relax":
            if random.uniform(0, 1) < 0.5:
                join_pose, idx_receptor, idx_ligand = self.define_relaxedbackbone(
                    ind, pose
                )
            else:
                join_pose, idx_receptor, idx_ligand = self.define_ensemble(ind, pose)
            return join_pose, idx_receptor, idx_ligand
        if bb_strategy == "popul_library":
            idx_receptor, idx_ligand = ind.idx_receptor, ind.idx_ligand
            pose_chainA = self.swap_operator.list_receptor[idx_receptor]
            pose_chainB = self.swap_operator.list_ligand[idx_ligand]
            join_pose = self.swap_operator.make_pose_with_chains(pose, pose_chainA, pose_chainB)
            return join_pose, idx_receptor, idx_ligand


        join_pose = pose
        idx_receptor, idx_ligand = ind.idx_receptor, idx_ligand
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
            apply_relax = self.relax_prob > random.uniform(0, 1)
            if not apply_relax:
                self.slide_into_contact.apply(join_pose)
                self.docking.apply(join_pose)
                after = self.energy_score(join_pose)
                # print("apply_dock_" + str(ind.idx) + " : " + str(after))
            else:
                self.slide_into_contact.apply(join_pose)
                self.docking.apply(join_pose)
                self.relax.apply(join_pose)
                remove_virtual_residues(join_pose)
                join_pose.fold_tree(self.native_fold_tree)
                after = self.energy_score(join_pose)
                if self.config.bb_strategy == "only_relax":
                    self.dock_pose = join_pose 
                else:
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
        if self.config.docking_type_option == "Unbound" and self.config.bb_strategy != "relax_best":
            return self.apply_unbound_docking(ind, local_search)
        else:
            return self.apply_bound_docking(ind, local_search)
