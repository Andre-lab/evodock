#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Author: Mads Jeppesen
@Date: 6/30/21
"""
import math
from copy import deepcopy
import time
import os
import copy
from itertools import islice
from pyrosetta.rosetta.protocols.moves import MonteCarlo # CycleMover #, TrialMover, SequenceMover, JumpOutMover, CycleMover
# from pyrosetta.rosetta.protocols.moves import MonteCarlo, CycleMover, TrialMover, SequenceMover, JumpOutMover, CycleMover
from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMover, MinMover, PackRotamersMover, TaskAwareMinMover
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.task_operations import RestrictToInterfaceVectorOperation
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, InitializeFromCommandline, IncludeCurrent,\
    NoRepackDisulfides, ReadResfile, AppendRotamerSet
from pyrosetta import Pose
from pyrosetta.rosetta.core.pack.rotamer_set import UnboundRotamersOperation
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.utility import vector1_unsigned_long
# from pyrosetta.rosetta.protocols.rigid import RigidBodyDofAdaptiveMover
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_jump_num, is_symmetric, jump_num_sym_dof
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.residue_selectors import TaskSelector
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, move_map_action
from pyrosetta.rosetta.protocols.relax import FastRelax
from cubicsym.actors.cubicboundary import CubicBoundary
from cubicsym.cubicsetup import CubicSetup
from cubicsym.dofspec import DofSpec
from cubicsym.kinematics import default_HF_dofs, default_HF_dofs_with_cubic_limits
from cloudcontactscore.cloudcontactscore import CloudContactScore
from cubicsym.actors.common_movers import TrialMover, JumpOutMover, SequenceMover, CycleMover
# from pyrosetta.rosetta.protocols.moves import MonteCarlo, CycleMover, TrialMover, SequenceMover, JumpOutMover, CycleMover
from cubicsym.actors.boundedminmover import BoundedMinMover
import numpy as np
from cubicsym.cubicmontecarlo import CubicMonteCarlo
from cloudcontactscore.cloudcontactscorecontainer import CloudContactScoreContainer
from pyrosetta.rosetta.basic.datacache import CacheableStringMap
from pyrosetta.rosetta.basic.datacache import CacheableStringFloatMap
from pyrosetta.rosetta.core.pose.datacache import CacheableDataType

def individual_is_within_bounds(config, fitnessfunc, ind):
    if config.syminfo is not None:
        pose = fitnessfunc.apply_genotype_to_pose(ind.genotype)
        return config.syminfo.cubicboundary.all_dofs_within_bounds(pose, raise_warning=True)
    else:
        pass

class SymShapeDock:

    def __init__(self, config, fafitnessfunc: "FAFitnessFunction", dock_attempts=10):
        self.config = config
        self.scfxn = fafitnessfunc.scfxn_rosetta
        self.dock_attempts = dock_attempts
        self.ccsc = self.config.syminfo.ccsc

        #self.idx_to_ccs, self.idx_to_cmc = None, None
        # # if running the same backbone all the time we only need 1 CloudContactScore/CubicMonteCarlo, else we need one for each backbone.
        # # if using multiple backbones we construct the CloudContactScore during apply
        # if config.flexbb:
        #     self.ccs, self.cmc = None, None
        #     self.idx_to_ccs, self.idx_to_cmc = {}, {}
        # else:
        #     self.set_ccs_and_cmc(fafitnessfunc.dock_pose)

        self.tf = self.create_taskfactory(fafitnessfunc.initial_pose)
        self.mc = MonteCarlo(self.scfxn, 0.8) # by NOT init pose here we have to set it later
        self.minmover = BoundedMinMover(self.config.syminfo.cubicboundary, self.scfxn)
        self.trial_minmover = TrialMover(self.minmover, self.mc)
        self.packer = self.create_packrotamers()
        self.trial_packer = TrialMover(self.packer, self.mc)
        self.trial_pack_and_mc_min = self.construct_pack_and_mc_min(fafitnessfunc.initial_pose)

    # def construct_ccs_and_cmc(self, pose):
    #     ccs = CloudContactScore(pose=pose, symdef=self.config.syminfo.input_symdef,
    #                                  use_atoms_beyond_CB=False, use_neighbour_ss=False)
    #     cmc = CubicMonteCarlo(scorefunction=ccs, cubicdofs=self.config.syminfo.cubicdofs)
    #     return ccs, cmc
    #
    # def set_ccs_and_cmc(self, pose):
    #     # get idx_subunit if present in the pose
    #     if pose.data().has(CacheableDataType.ARBITRARY_STRING_DATA):
    #         assert self.config.flexbb is True
    #         stringmap = pose.data().get_ptr(CacheableDataType.ARBITRARY_STRING_DATA)
    #         idx_subunit = stringmap.map()["idx_subunit"] # Must be defined in this case!
    #     else:
    #         assert self.config.flexbb is False
    #         idx_subunit = None
    #     # if no idx_subunit is defined, then construct from new
    #     if idx_subunit is None:
    #         self.ccs, self.cmc = self.construct_ccs_and_cmc(pose)
    #     else:
    #         # if idx_subunit is defined, then check if it has already been created.
    #         # if it has, then return them, if not, construct them and save them.
    #         if idx_subunit in self.idx_to_ccs:
    #             self.ccs = self.idx_to_ccs[idx_subunit]
    #             self.cmc = self.idx_to_cmc[idx_subunit]
    #         else:
    #             print(f"Constructing ccs and cmc for {idx_subunit}")
    #             self.ccs, self.cmc = self.construct_ccs_and_cmc(pose)
    #             self.idx_to_ccs[idx_subunit] = self.ccs
    #             self.idx_to_cmc[idx_subunit] = self.cmc

    @staticmethod
    def create_taskfactory(pose, cb_dist_cutoff: int = 10.0, nearby_atom_cutoff: int = 5.5,
                           vector_angle_cutoff: int = 75.0, vector_dist_cutoff: int = 9.0,
                           include_all_water: bool = False) -> TaskFactory:
        """Creates a taskfactory. Attempting to imitate protocols.docking.DockTaskFactory

        :param pose: Pose to create taskfactory from.
        :param cb_dist_cutoff: Option for RestrictToInterfaceVector.
        :param nearby_atom_cutoff: Option for RestrictToInterfaceVector.
        :param vector_angle_cutoff: Option for RestrictToInterfaceVector.
        :param vector_dist_cutoff: Option for RestrictToInterfaceVector.
        :param include_all_water: Option for RestrictToInterfaceVector.
        :return: TaskFactory
        """
        taskfactory = TaskFactory()
        # I am using the Restricttointerfacevector vs RestrictToInterface, but both I believe would work.
        pose.num_chains()
        task1 = RestrictToInterfaceVectorOperation()
        # assume chain A is the upper_chain and the lower_chain are the others
        task1.lower_chain(1)
        chain_num_upper = vector1_unsigned_long()
        for chain in range(2, pose.num_chains() + 1):
            chain_num_upper.append(chain)
        task1.upper_chain(chain_num_upper)
        # set other setters with default values
        task1.CB_dist_cutoff(cb_dist_cutoff)
        task1.nearby_atom_cutoff(nearby_atom_cutoff)
        task1.vector_angle_cutoff(vector_angle_cutoff)
        task1.vector_dist_cutoff(vector_dist_cutoff)
        task1.include_all_water(include_all_water)
        taskfactory.push_back(task1)
        # -- create the once already existing in the original DockMCMProtocol -- #
        taskfactory.push_back(RestrictToRepacking())
        taskfactory.push_back(InitializeFromCommandline())
        taskfactory.push_back(IncludeCurrent())
        taskfactory.push_back(NoRepackDisulfides())
        taskfactory.push_back(ReadResfile())
        unboundrot = UnboundRotamersOperation()
        unboundrot.initialize_from_command_line()
        taskfactory.push_back(AppendRotamerSet(unboundrot))
        return taskfactory

    def create_packrotamers(self):
        """Returns a PackRotamersMover"""
        packer = PackRotamersMover(self.scfxn)
        packer.task_factory(self.tf)
        # packer = RotamerTrialsMover(self.scfxn, self.tf)
        return packer

    def construct_pack_and_mc_min(self, pose):
        """Returns a JumpOutMover that applies a PackRotamersMover and energy dependingly a BoundedMinMover."""
        #multiplier = self.config.syminfo.cubicboundary.symmetrysetup.cubic_energy_multiplier_from_pose(pose)
        minimization_threshold = 15.0 #* multiplier
        rb_mover_min = JumpOutMover(self.packer, self.minmover, self.scfxn, minimization_threshold)
        return TrialMover(rb_mover_min, self.mc)

    def set_constraints(self, pose):
        self.config.syminfo.cubicboundary.set_constraints(pose)
        assert self.config.syminfo.cubicboundary.constrains_are_set_in_score(
            self.scfxn), "In order for constraints to work the weights need to be set in the scorefunction"

    # todo: instead of random search use mc/greedyalgo to go towards the least clashes?
    # todo: check if the good low score of the old behavior is because it goes out of symmetry?
    # todo: check if this is how you should use the rb_mover as below or should you call reset?
    # todo: check if you should use the multiplier (x60)

    def reset(self, pose):
        """Resets certain objects. Important before any moves are done!"""
        self.mc.reset(pose) # reset montecarlo
        self.ccsc.cmc.reset(pose) # reset cubicmontecarlo
        # self.rb_mover.reset_pose(pose)  # records the initial placement of the pose. Used to perturb the dofs.
        # self.rb_mover.reset()  # resets the params to the initial values
        self.rb_mover = self.config.syminfo.cubicboundary.construct_rigidbody_mover(pose)

    # fixme: I do not understand if self.mc.mc_accepted() == 1 fully and therefore will revert to use self.scfxn directly
    def get_current_pose_score_quick(self, pose):
        """Gets the current score of the pose. Bypasses calling the score function directly as the current score of the pose
        is stored in the MonteCarlo object. This should be faster."""
        return self.scfxn.score(pose)
        # if self.mc.last_accept():
        # if self.mc.mc_accepted() == 1:
        #     score = self.mc.last_score()
        # else: # self.mc.mc_accepted() == 2:
        #     score = self.mc.lowest_score()
        # # elif self.mc.mc_accepted() == 1:
        # # fixme Delete this assert!
        # assert math.isclose(score, self.scfxn.score(pose)), "energy does not match!"
        # return score

    def apply(self, pose):
        # set the correct ccs and cmc
        self.ccsc.set_ccs_and_cmc(pose)
        self.set_constraints(pose)
        self.reset(pose)
        # -- 1: initial pack --
        # calc score for later delta
        pre_score = self.get_current_pose_score_quick(pose)
        self.trial_packer.apply(pose)
        # -- 2: initial min --
        if self.get_current_pose_score_quick(pose) - pre_score < 15:
            self.trial_minmover.apply(pose)
        # calc score for later delta
        pre_score = self.get_current_pose_score_quick(pose)
        # -- 3: docking cycle
        self.inner_protocol(pose)
        # -4: final packing
        self.trial_packer.apply(pose)
        # 5: final min (if delta E is below 15)
        if self.get_current_pose_score_quick(pose) - pre_score < 15:
            self.trial_minmover.apply(pose)
        # 6 return the best score pose
        self.mc.recover_low(pose)

    def inner_protocol(self, pose):
        for it in range(self.dock_attempts):
            self.rb_mover.apply(pose)
            self.ccsc.cmc.apply(pose)
        self.ccsc.cmc.recover_lowest_scored_pose(pose)

    # def inner_protocol(self, pose):
    #     got_0_clashes, n_applies = False, 0
    #     # Try docking attempts
    #     best_positions = self.config.syminfo.cubicboundary.cubicdofs.get_positions_as_list(pose)
    #     for _ in range(self.dock_attempts):
    #         if n_applies == self.max_applies:
    #             break
    #         # apply dock move
    #         self.rb_mover.apply(pose)
    #         # record clashes
    #         n_clashes = self.ccs.number_of_clashes(pose)
    #         self.clashes_2_dofs[n_clashes] = self.config.syminfo.cubicboundary.cubicdofs.get_positions_as_list(pose)
    #         # if dock move gave 0 clashes, then apply packing/minimization. Else go back to the original position
    #         if n_clashes == 0:
    #             self.trial_pack_and_mc_min.apply(pose)
    #             got_0_clashes = True
    #             n_applies += 1
    #             best_positions = self.config.syminfo.cubicboundary.cubicdofs.get_positions_as_list(pose)
    #         else:
    #             self.config.syminfo.cubicboundary.cubicdofs.transfer_dofs_to_pose(pose, *best_positions)
    #     if not got_0_clashes:
    #         # find the solution with the lowest amount of clashes:
    #         clashes, best_sol = min(self.clashes_2_dofs.items(), key=lambda x: x[0])
    #         self.config.syminfo.cubicboundary.cubicdofs.transfer_dofs_to_pose(pose, *best_sol)
    #         self.trial_pack_and_mc_min.apply(pose)

class SymInfo:
    """Class for storing and accessing symmetrical information"""

    def __init__(self):
        """Initialize object"""
        self.input_symdef = None
        self.native_symdef = None
        self.jumps_str = None
        self.dofs_str = None
        self.bounds = None
        self.genotype_size = None
        self.initialize_rigid_body_dofs = None
        self.normalize_trans = None
        self.cubic_limits = None
        self.normalize_trans_map = None
        self.cubicboundary = None
        self.org_normalize_trans_map = None
        self.jumps_int = None
        self.dofs_int = None
        self.initial_placement = None
        self.ccsc = None
        self.bound_penalty = None
        self.native_symmetric_input = None
        # if pose:
        #     self.store_info_from_pose(pose)
        # if config:
        #     self.store_info_from_config(config)

    def get_position_info(self, pose: Pose) -> list:
        return self.dof_spec.get_positions_as_list(pose)

    # def store_info_from_config(self, config):
    #     assert config.has_option("Inputs", "symdef_file")
    #     # self.input_symdef = os.getcwd() + config.get("Inputs", "symdef_file")
    #     #self.native_symdef = os.getcwd() + config.get("Native", "symdef_file")
    #     if config.has_option("Bounds", "symdofs"):
    #         self.jumps_str = [i.split(":")[0] for i in config.get("Bounds", "symdofs").split(",")]
    #         self.dofs_str = [i.split(":")[1:] for i in config.get("Bounds", "symdofs").split(",")]
    #         bounds = iter([i for i in config.get("Bounds", "bounds").split(",")])
    #         # for nicer mapping we have to map the bounds the same ways as the dofs
    #         self.bounds = [list(islice(bounds, l)) for l in [len(i) for i in self.dofs_str]]
    #         self.genotype_size = config.get("Bounds", "symdofs").count(":")
    #         assert len(self.jumps_str) == len(self.dofs_str)
    #         assert len(self.jumps_str) == len(self.dofs_str)
    #     else:  # apply defaults
    #         raise NotImplementedError
    #     if config.has_option("Bounds", "init"):
    #         init_bounds = iter([i for i in config.get("Bounds", "init").split(",")])
    #         init_bounds = [list(islice(init_bounds, l)) for l in [len(i) for i in self.dofs_str]]
    #         self.init_bounds = [(-l, l) for l in [float(i[0]) / float(b[0]) for b, i in zip(self.bounds, init_bounds)]]
    #     else:
    #         self.init_bounds = None
    #     self.normalize_trans = [int(i) for i in config.get("Bounds", "normalize_trans", fallback="2000,1000").split(",")]
    #     self.bound_penalty = config.getfloat("Bounds", "bound_penalty", fallback=1)


    def store_info_from_pose(self, pose):
        """Store symmetry information from pose."""
        self.dof_spec = DofSpec(pose)
        self.dof_spec.set_symmetrical_bounds(self.bounds)
        self.genotype_size = self.dof_spec.dofsize
        assert self.bound_penalty is not None
        self.cubicboundary = CubicBoundary(CubicSetup(symdef=self.input_symdef), pose_at_initial_position=pose, dof_spec=self.dof_spec, sd=self.bound_penalty)
        self.initial_placement = self.dof_spec.get_positions_as_list(pose)
        self._map_normalize_trans_to_jumpdofs()

    def _map_normalize_trans_to_jumpdofs(self) -> None:
        if self.normalize_trans:
            trans = deepcopy(self.normalize_trans)
            self.normalize_trans_map = []
            for jump_int, dof_int in self.dof_spec.get_translational_dofs_int():
                transmag = trans.pop(0)
                self.normalize_trans_map.append([jump_int, dof_int, transmag])
            self.org_normalize_trans_map = copy.deepcopy(self.normalize_trans_map)