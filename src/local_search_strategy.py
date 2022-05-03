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
from pyrosetta.rosetta.protocols.symmetry import FaSymDockingSlideTogether
from src.symmetry import SymDockMCMProtocol, SymDockingSlideIntoContactWrapper, SequentialSymmetrySliderWrapper
from pyrosetta.rosetta.core.conformation.symmetry import SlideCriteriaType
from pyrosetta.rosetta.protocols.symmetry import SymmetrySlider
from pyrosetta.rosetta.protocols.symmetry import SequentialSymmetrySlider
from pyrosetta.rosetta.protocols.moves import NullMover
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from src.utils import IP_ADDRESS


# # <<<<<<< HEAD
#     # Options:
#     # None: only score and return the poses
#     # only_slide: just slide_into_contact
#     # mcm_rosetta: mcm protocol mover (high res) from rosetta (2 cycles)
#     def __init__(self, scfxn, packer_option="default_combination", slide=True, show_local_search=False,
#         pymol_history=False):
#         self.packer_option = packer_option
#         self.scfxn = scfxn
#         self.local_logger = logging.getLogger("evodock.local")
#         self.local_logger.setLevel(logging.INFO)
#         self.slide = slide
#         self.show_local_search = show_local_search
#         self.pymol_history = pymol_history
#
#         if is_symmetric(scfxn.dock_pose):
#             self.slide_into_contact = SequentialSymmetrySlider(scfxn.dock_pose, SlideCriteriaType(1))
#
#             # This will only slide on the first pose!!
#             # FA_REP_SCORE = SlideCriteriaType(2)
#             # self.slide_into_contact = SequentialSymmetrySlider(scfxn.dock_pose, FA_REP_SCORE)
#
#             # todo: use the highresolution alternative below or delete it
#             # dofs = scfxn.dock_pose.conformation().Symmetry_Info().get_dofs()
#             # self.slide_into_contact = FaSymDockingSlideTogether(dofs)
#         else:
#             self.slide_into_contact = DockingSlideIntoContact(1)
#         if packer_option == "mcm_rosetta":
#             if is_symmetric(scfxn.dock_pose):
#                 self.docking = SymDockMCMProtocol(scfxn.dock_pose)
#             else:
#                 mcm_docking = DockMCMProtocol()
#                 mcm_docking.set_native_pose(scfxn.dock_pose)
#                 mcm_docking.set_scorefxn(scfxn.scfxn_rosetta)
#                 mcm_docking.set_rt_min(False)
#                 mcm_docking.set_sc_min(False)
#                 mock_pose = Pose()
#                 mock_pose.assign(scfxn.dock_pose)
#                 mcm_docking.apply(mock_pose)
#                 self.docking = mcm_docking
#                 # DEBUG Is this a hack to skip the create_and_attach_task_factory on each apply to save time??
#                 # A reason for potentially deleting this is that it is quite confusing. The taskfactory is already set when
#                 # calling apply above (see protocols.docking.DockMCMProtocol.cc:179). It is the default one in this case.
#                 # Then the taskfacory is set again using the default one, in a way that is intented for a new tasks using
#                 # task_factory(). Although this seems to actually surpass the creation of  the task_factory on each apply
#                 # so it could be considered smart.
#                 self.docking.set_task_factory(mcm_docking.task_factory())
#                 self.docking.set_ignore_default_task(True)
# # =======

class LocalSearchStrategy:
    def __init__(self, config, scfxn, dock_pose):
        # self.pymover = PyMOLMover(address="10.8.0.22", port=65000, max_packet_size=1400)
        self.config = config
        self.scfxn = scfxn
        self.dock_pose = dock_pose
        self.packer_option = config.local_search_option
        self.native_fold_tree = scfxn.dock_pose.fold_tree()
        # todo for symmetry!
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
            if is_symmetric(scfxn.dock_pose):
                self.docking = SymDockMCMProtocol(scfxn.dock_pose)
            else:
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
        # slide option
        if self.config.slide:
            if is_symmetric(scfxn.dock_pose):
                self.slide_into_contact = SequentialSymmetrySlider(scfxn.dock_pose, SlideCriteriaType(1))
                # This will only slide on the first pose!!
                # FA_REP_SCORE = SlideCriteriaType(2)
                # self.slide_into_contact = SequentialSymmetrySlider(scfxn.dock_pose, FA_REP_SCORE)
                # todo: use the highresolution alternative below or delete it
                # dofs = scfxn.dock_pose.conformation().Symmetry_Info().get_dofs()
                # self.slide_into_contact = FaSymDockingSlideTogether(dofs)
            else:
                self.slide_into_contact = FaDockingSlideIntoContact(dock_pose.num_jump())
        else:
            self.slide_into_contact = NullMover()
        if self.config.docking_type_option == "Flexbb":
            self.swap_operator = FlexbbSwapOperator(config, scfxn, None)

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def apply_bb_strategy(self, ind, pose):
        idx_receptor, idx_ligand, idx_subunit = ind.idx_receptor, ind.idx_ligand, ind.idx_subunit
        if self.config.syminfo:
            join_pose = self.swap_operator.list_subunits[ind.idx_subunit]
            for jump in self.config.syminfo.get("jumps_int"):
                join_pose.set_jump(jump, pose.jump(jump))
        else:
            pose_chainA = self.swap_operator.list_receptor[idx_receptor]
            pose_chainB = self.swap_operator.list_ligand[idx_ligand]
            join_pose = self.swap_operator.make_pose_with_chains(
                pose, pose_chainA, pose_chainB
            )
        return join_pose, idx_receptor, idx_ligand, idx_subunit

        # join_pose = pose
        # idx_receptor, idx_ligand = ind.idx_receptor, idx_ligand
        # return join_pose, idx_receptor, idx_ligand

    def apply_bound_docking(self, ind, local_search=True):
        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
        before = self.energy_score(pose)
        if local_search and self.packer_option != "None":
            if self.config.show_local_search:
                pose.pdb_info().name(f"IND{ind.idx}")
                self.config.pmm.apply(pose)
            self.slide_into_contact.apply(pose)
            self.docking.apply(pose)
            after = self.energy_score(pose)
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
        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
        # pose.pdb_info().name("apply_gen_" + str(ind.idx))
        # self.pymover.apply(pose)
        join_pose, idx_receptor, idx_ligand, idx_subunit = self.apply_bb_strategy(ind, pose)
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
            "idx_subunit": idx_subunit,
        }
        return return_data

    def apply(self, ind, local_search=True):
        if self.config.docking_type_option == "Flexbb":
            return self.apply_unbound_docking(ind, local_search)
        else:
            return self.apply_bound_docking(ind, local_search)
