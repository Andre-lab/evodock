#!/usr/bin/env python
# coding: utf-8
from src.flexbb_swap_operator import FlexbbSwapOperator
# from pyrosetta.rosetta.protocols.moves import PyMOLMover
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
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pyrosetta.rosetta.protocols.moves import NullMover
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta import pose_from_file
from cubicsym.actors.cubicsymmetryslider import CubicSymmetrySlider
from src.symmetry import SymDockMCMProtocol, DockNRelaxProtocol, SymShapeDock
from pyrosetta.rosetta.core.kinematics import Jump
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.core.pose.symmetry import extract_asymmetric_unit

from symmetryhandler.reference_kinematics import set_all_dofs_to_zero
from cloudcontactscore.cloudcontactscorecontainer import CloudContactScoreContainer
from cubicsym.utilities import add_base_to_pose, add_id_to_pose_w_base


class LocalSearchStrategy:
    def __init__(self, config, scfxn, dock_pose):
        # self.pymover = PyMOLMover(address="10.8.0.22", port=65000, max_packet_size=1400)
        self.config = config
        self.scfxn = scfxn
        self.dock_pose = dock_pose
        self.packer_option = config.local_search_option
        self.native_fold_tree = scfxn.dock_pose.fold_tree()
        if self.config.syminfo:
            self.symmetrymover = SetupForSymmetryMover(self.config.syminfo.input_symdef)
        else:
            self.symmetrymover = None
        # docking option
        if is_symmetric(scfxn.dock_pose):
            # construct a CloudContactScoreContainer in the SymInfo object
            self.config.syminfo.ccsc = CloudContactScoreContainer(self.dock_pose, self.config.syminfo.input_symdef,
                                                                  low_memory= self.config.low_memory_mode)
            if self.packer_option == "symshapedock":
                self.docking = SymShapeDock(config, scfxn)
            else:
                raise NotImplementedError(f"Only 'symshapedock' is a valid local_search options")
        else:
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
            elif self.packer_option == "mcm_rosetta":
                mcm_docking = DockMCMProtocol()
                mcm_docking.set_first_cycle(self.config.num_first_cycle)
                mcm_docking.set_second_cycle(self.config.num_second_cycle)
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
            else:
                raise NotImplementedError(f"Only 'sidechains', 'mcm_rosetta' are valid local_search options")
        # slide option
        if self.config.slide:
            if is_symmetric(scfxn.dock_pose):
                self.slide_into_contact = CubicSymmetrySlider(dock_pose, self.config.syminfo.input_symdef, ccsc=self.config.syminfo.ccsc,
                                                              trans_mag=self.config.slide_trans_mag, pymolmover=self.config.pmm,
                                                              max_slide_attempts=self.config.max_slide_attempts,
                                                              cubicboundary=self.config.syminfo.cubicboundary)
            else:
                self.slide_into_contact = FaDockingSlideIntoContact(dock_pose.num_jump())
        else:
            self.slide_into_contact = NullMover()
        if self.config.flexbb:
            self.swap_operator = FlexbbSwapOperator(config, scfxn, self)

    def extract_single_chain_and_reset_dof(self, pose):
        """On a copy of the pose, resets all Jumps in the (symmetrical) pose to unity and return the main chain."""
        assert is_symmetric(pose)
        pose_reset = pose.clone()
        # for jump in self.config.syminfo.jumps_int:
        #     pose_reset.set_jump(jump, Jump()) # set with emtpy jumps
        set_all_dofs_to_zero(pose_reset)
        out_pose = Pose()
        extract_asymmetric_unit(pose_reset,out_pose, False)
        return out_pose

    def superimpose_tobe_symmetrical_pose(self, pose, join_pose):
        # reset all the dofs for the pose to be aligned to
        pose_reset = self.extract_single_chain_and_reset_dof(pose)
        # align it
        calpha_superimpose_pose(join_pose, pose_reset)

    def apply_bb_strategy(self, ind, pose):
        idx_receptor, idx_ligand, idx_subunit = ind.idx_receptor, ind.idx_ligand, ind.idx_subunit
        if self.config.syminfo:
            join_pose = self.swap_operator.list_subunits[ind.idx_subunit]
            if self.config.low_memory_mode:
                join_pose = pose_from_file(join_pose)
            else:
                # cloning is to make sure we don't symmetrize the stored pose in self.list_subunits because this can swamp the memory
                assert not is_symmetric(join_pose), "join_pose should not be stored symmetrically as it hocks the memory!"
                join_pose.clone()
            # we have to align join_pose and pose onto one another
            self.superimpose_tobe_symmetrical_pose(pose, join_pose)
            # symmetrize the join_pose
            self.symmetrymover.apply(join_pose)
            # calpha_superimpose_pose
            for jump in self.config.syminfo.dof_spec.jump_int:
                join_pose.set_jump(jump, pose.jump(jump))
        else:
            pose_chainA = self.swap_operator.list_receptor[idx_receptor]
            pose_chainB = self.swap_operator.list_ligand[idx_ligand]
            join_pose = self.swap_operator.make_pose_with_chains(
                pose, pose_chainA, pose_chainB
            )
        return join_pose, idx_receptor, idx_ligand, idx_subunit

    def apply_bound_docking(self, ind, local_search=True):
        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
        before = self.scfxn.score(pose)
        if local_search and self.packer_option != "None":
            if self.config.show_local_search:
                pose.pdb_info().name(f"IND{ind.idx}")
                self.config.pmm.apply(pose)
            if self.config.syminfo is not None:
                add_base_to_pose(pose)
            self.slide_into_contact.apply(pose)
            self.docking.apply(pose)
            after = self.scfxn.score(pose)
            if self.config.show_local_search:
                pose.pdb_info().name(f"IND{ind.idx}")
                self.config.pmm.apply(pose)
        else:
            after = before

        return_data = {
            "pose": pose,
            "before": before,
            "after": after,
            "idx_ligand": ind.idx_ligand,
            "idx_receptor": ind.idx_receptor,
            "idx_subunit": ind.idx_subunit
        }
        return return_data


    def apply_unbound_docking(self, ind, local_search=True):
        # Generate the current pose with the genotype, but it does not have the correct backbone as it only uses the ind.genotype
        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
        self.config.visualize_pose(pose, ind.idx)
        # Generate the correct backbone as this, internally uses the ind.idx_<>
        join_pose, idx_receptor, idx_ligand, idx_subunit = self.apply_bb_strategy(ind, pose)
        before = self.scfxn.score(join_pose)
        if local_search and self.packer_option != "None":
            if self.config.syminfo is not None:
                add_id_to_pose_w_base(join_pose, idx_subunit)
            self.slide_into_contact.apply(join_pose)
            self.docking.apply(join_pose)
            after = self.scfxn.score(join_pose)
        else:
            after = before
        return_data = {
            "pose": join_pose,
            "before": before,
            "after": after,
            "idx_ligand": idx_ligand,
            "idx_receptor": idx_receptor,
            "idx_subunit": idx_subunit,
        }
        return return_data

    def apply(self, ind, local_search=True):
        if self.config.flexbb:
            return self.apply_unbound_docking(ind, local_search)
        else:
            return self.apply_bound_docking(ind, local_search)
