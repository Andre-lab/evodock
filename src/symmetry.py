#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Symmetry classes: SymDockMCMProtocol and SymDockMCMCycle

To my knowledge there's no equivalent symmetrical version of DockMCMProtocol in Rosetta. The 2 classes describe here is an
attempt at reproducing the same behavior of DockMCMProtocol in a symmetrical form and as it is used in EvoDOCK from the Rosetta C++ code.
The following problems using the DockMCMProtocol for symmetry are:

1. The use of the RigidBodyPerturbMover.
   - It does not consider symmetry and so will mess it up if used.
2. The use of DockingNoRepack.
   - It seems like the  mover also doesnt work for symmetry. It defines only a single jump.

TODO: If one want to implement sc or bb minimzation during the protocol the protocol has to be implemented. Some of the code is already there
TODO: if one want to constrain the RigidBodyDofAdaptive in the future - one have to call reset() and put initialization in the constructor

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
from pyrosetta.rosetta.protocols.rigid import RigidBodyDofAdaptiveMover
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
    pose = fitnessfunc.apply_genotype_to_pose(ind.genotype)
    return config.syminfo.cubicboundary.all_dofs_within_bounds(pose, raise_warning=True)

class CycleWClashMover:
    """Class with same function as the CycleMover in Rosetta but skips the following steps in the cycle if clashes are detected."""

    def __init__(self, iterations, ccs: CloudContactScore):
        """

        :param iterations: Number of iterations to use
        :param ccs: CloudContactScore object ot use.
        :return: None
        """
        self.iterations = iterations
        self.ccs = ccs
        self.movers = []
        self.clashcheck = []

    def add_clash_check(self):
        """Add a clash check"""
        self.movers.append(None)
        self.clashcheck.append(True)

    def add_move(self, mover):
        """Add a mover"""
        self.movers.append(mover)
        self.clashcheck.append(False)

    def apply(self, pose):
        """Apply the cycle on the pose"""
        for it in range(self.iterations):
            for clashcheck, mover in zip(self.clashcheck, self.movers):
                if not clashcheck or (clashcheck and self.ccs.number_of_clashes(pose) == 0):
                    mover.apply(pose)
                else:
                    break

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

        self.tf = self.create_taskfactory(fafitnessfunc.dock_pose)
        self.mc = MonteCarlo(self.scfxn, 0.8) # by NOT init pose here we have to set it later
        self.minmover = BoundedMinMover(self.config.syminfo.cubicboundary, self.scfxn)
        self.trial_minmover = TrialMover(self.minmover, self.mc)
        self.packer = self.create_packrotamers()
        self.trial_packer = TrialMover(self.packer, self.mc)
        self.trial_pack_and_mc_min = self.construct_pack_and_mc_min(fafitnessfunc.dock_pose)

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

class SymDockMCMCycle:
    """The symmetric version of DockMCMCycler class in C++ with extra functionality. The extra functionality comes from the posibility
    of specifying a dofspecification that fine-tunes seperate symmetric dofs. See the comments in setup_protocol for a
    description of the algorithm."""

    def __init__(self, pose, config, scorefxn, dofspecification: dict, repack_period=8):
        """Initialization"""
        self.scorefxn = scorefxn
        self.set_default(pose, repack_period, dofspecification, minimization_option=config.minimization_option,
                         min_option_bb=config.min_option_bb, min_option_sc=config.min_option_sc, cartesian=config.cartesian)
        self.dock_mcm_cycle = None
        self.config = config

    def reset_cycle_index(self):
        """Resets the internal cycle mover to step 0 (next_move_ == 0 in C++)"""
        self.dock_mcm_cycle.reset_cycle_index()

    def apply(self, pose):
        """Applies docking. See the comments in setup_protocol for a description of the algorithm"""
        # initialize the protocol if it hasn't been set
        if not self.dock_mcm_cycle:
            self.init_mc(pose)
        # Apply the protocol
        self.dock_mcm_cycle.apply(pose)

    def init_mc(self, pose):
        """Initializes the MonteCarlo object and sets up the dock_mcm_cycle objet to be used in the 'apply' function."""
        # in the c++ code it is written that this is to ensure that consistency if
        # the mover is shared across multiple other movers
        if self.mc.last_accepted_pose().empty():
            self.mc.reset(pose)
        self.setup_protocol(pose)

    def set_default(self, pose, repack_period, dof_specification, sc_min=False, rt_min=False, bb_min_res:list=None, sc_min_res:list=False, minimization_option="jumps",
                    min_option_bb=False, min_option_sc=False, cartesian=False):
        """Sets up the default parameters used when creating the dock_mcm_cycle object in the 'setup_protocol' function.

        :param pose:
        :param repack_period: basically the cycles of the docking protocol. Default is 8.
        :param dof_specification: specifying the jumpdof parameters
        :param sc_min: Have a SidechainMinMover at the end of the docking protocol cycle.
        :param rt_min: Have a RotamertrialsMinMover at the end of the docking protocol cycle.
        :param bb_min_res: Turns on the given bb residue positions in the movemap that is always used by the internal MinMover.
        :param sc_min_res: Turns on the given sc residue positions in the movemap that is always used by the internal MinMover.
        :return: None
        """
        self.repack_period = repack_period
        self.dofspecification = dof_specification
        self.sc_min = sc_min
        self.rt_min = rt_min
        self.cartesian = cartesian
        if self.cartesian:
            assert self.scorefxn.get_name() == "ref2015", "Code doesnt know how to create a cartesian version of the scorefxn - only ref2015 is supported."
            self.scorefxn_cart = ScoreFunctionFactory.create_score_function("ref2015_cart")
            self.config.syminfo.cubicboundary.turn_on_constraint_for_score(self.scorefxn_cart)
        self.min_tolerance = 0.01
        self.min_type = "lbfgs_armijo_nonmonotone"
        self.nb_list = True

        # DockMCM mover only initializes an empty taskfatory so one have to set it manually here
        self.tf = self.create_taskfactory(pose)
        self.mc = MonteCarlo(pose, self.scorefxn, 0.8)
        self.repack_period = repack_period
        self.movemapfactory, self.movemap, self.min_option_bb, self.min_option_sc = None, None, min_option_bb, min_option_sc
        assert minimization_option in ("jumps", "interface"), minimization_option
        if minimization_option == "interface":
            self.movemapfactory = self.create_interface_mf()
        elif minimization_option == "jumps": # classic DockMCMMover where only the jumps are minimized
            self.movemap = MoveMap()
            self.movemap.set_chi(False)
            self.movemap.set_bb(False)
            # Jumps are always minimized
            for jump_name, _ in self.dofspecification.items():
                jumpid = sym_dof_jump_num(pose, jump_name)
                self.movemap.set_jump(jumpid, True)
            # Choose to minimize parts of the bb and sc
            if bb_min_res:
                for pos in bb_min_res:
                    self.movemap.set_bb(pos, True)
            if sc_min_res:
                for pos in sc_min_res:
                    self.movemap.set_bb(pos, True)
        assert self.movemapfactory or self.movemap and not (self.movemapfactory != None and self.movemap != None), \
            "either a movemapfactory or movemap has to be parsed and not both!"


    # if bb_min_res and sc_min_res:
    #     exit("bb_min_res and sc_min_res cannot be used when a minimization_option is being used.")
    # if minimization_option == "interface":

    def setup_protocol(self, pose):
        """Sets up the dock_mcm_cycle object."""
        # Contruction of rb_mover_min_trial (The first move in the sequence)
        # 1: (1a. rb dock move -> 1b. rotamer trial) -> 2: if E < 15.0 apply minimization.
        # Accept the entire move 1+2 with the MC criterion
        self.rb_mover = self.config.syminfo.cubicboundary.construct_rigidbody_mover(pose)
        rottrial = RotamerTrialsMover(self.scorefxn, self.tf)
        rb_pack_min = SequenceMover()
        rb_pack_min.add_mover(self.rb_mover)
        rb_pack_min.add_mover(rottrial)
        minimization_threshold = 15.0
        min_mover = self.construct_minmover()
        rb_mover_min = JumpOutMover(rb_pack_min, min_mover, self.scorefxn, minimization_threshold)
        rb_mover_min_trial = TrialMover(rb_mover_min, self.mc)

        # Contruction of repack_step (The Second move in the sequence)
        # 1: Exactly the same as rb_mover_min_trial (1: (1a. rb dock move -> 1b. rotamer trial) -> 2: if E < 15.0 apply minimization.)
        # 2: Packrotamers
        # Accept the entire move 1+2 with the MC criterion
        repack_step = SequenceMover()
        repack_step.add_mover(rb_mover_min_trial)
        packer = PackRotamersMover(self.scorefxn)
        packer.task_factory(self.tf)
        pack_interface_and_move_loops_trial = TrialMover(packer, self.mc)
        repack_step.add_mover(pack_interface_and_move_loops_trial)
        if self.rt_min:
            pass
        if self.sc_min:
            pass
        dock_mcm_cycle = CycleMover()
        for i in range(self.repack_period):
            dock_mcm_cycle.add_mover(rb_mover_min_trial)
            dock_mcm_cycle.add_mover(repack_step)

        self.dock_mcm_cycle = dock_mcm_cycle

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

    def construct_minmover(self):
        """Wrapper for the initialization of a MinMover."""
        minmover = MinMover()
        minmover.cartesian(self.cartesian)
        if self.cartesian:
            minmover.score_function(self.scorefxn_cart)
        else:
            minmover.score_function(self.scorefxn)
        minmover.tolerance(self.min_tolerance)
        minmover.min_type(self.min_type)
        minmover.nb_list(self.nb_list)
        # if minimize_option == "interface":
        if self.movemapfactory:
            minmover.movemap_factory(self.movemapfactory)
        else: # regular minimization as DockMCMMover
            minmover.movemap(self.movemap)
        return minmover

    def create_interface_mf(self):
        # construct mf and set all jumps to true
        mf = MoveMapFactory()
        mf.all_jumps(True)
        # create residueselector that just selects the interface
        residueselector = TaskSelector()
        residueselector.set_select_designable(True)
        residueselector.set_select_packable(True)
        residueselector.set_select_fixed(False)
        residueselector.set_task_factory(self.tf)
        if self.min_option_bb:
            mf.add_bb_action(move_map_action.mm_enable, residueselector)
        if self.min_option_sc:
            mf.add_chi_action(move_map_action.mm_enable, residueselector)
        # just make sure no bondangles/lengths are minimized
        mf.all_bondangles(False)
        mf.all_bondlengths(False)
        return mf

class SymDockMCMProtocol:
    """The symmetric version of DockMCMProtocol class in C++ with extra functionality."""

    def __init__(self, config, fafitnessfunc: "FAFitnessFunction"):
        """Initializes the mover.

        :param pose: Pose to dock.
        :param num_of_first_cycle: Number of iterations for the first cycle
        :param num_of_second_cycle: Number of iterations for the second cycle
        """
        self.num_of_first_cycle = config.num_first_cycle
        self.num_of_second_cycle = config.num_second_cycle
        self.scorefxn = fafitnessfunc.scfxn_rosetta
        assert config.syminfo.cubicboundary.constrains_are_set_in_score(self.scorefxn), "In order for constraints to work these need to be set in the scorefunction"
        # these are form DockingHighRes - they are currently not used in the protocol
        self.sc_min = False
        self.rt_min = False
        self.config = config
        # todo: has to be changed for flexbb
        assert config.flexbb is False
        self.ccs = CloudContactScore(pose=fafitnessfunc.dock_pose, symdef=config.syminfo.input_symdef)

        # todo: delete this:
        self.initial_pack_time = []
        self.first_min_time = []
        self.dock_mcm_time = []
        self.final_min_time = []

    def apply(self, pose):
        """Applies the mover."""
        assert is_symmetric(pose), "Pose is not symmetric!"

        if self.config.old_behavior_w_cont:
            self.config.syminfo.cubicboundary.set_constraints(pose)
            assert self.config.syminfo.cubicboundary.constrains_are_set_in_score(self.scorefxn), "In order for constraints to work these need to be set in the scorefunction"

        dock_mcm = SymDockMCMCycle(pose, self.config, self.scorefxn, self.config.syminfo.cubicboundary.dofspecification)

        import time
        start = time.time()
        initial_pack = PackRotamersMover()
        initial_pack.score_function(self.scorefxn)
        initial_pack.task_factory(dock_mcm.tf)
        initial_pack_trial = TrialMover(initial_pack, dock_mcm.mc)
        initial_pack_trial.apply(pose)
        self.initial_pack_time.append(time.time() - start)
        print("initial_pack = time:", np.array(self.initial_pack_time).mean())

        # 2. Minimization, accept with MC
        # Called DockMinMover in C++, but this just seems to be a minmover in a trialmover with the exact same minmover settings as dock_mcm
        start = time.time()
        if self.config.old_behavior:
            minmover = dock_mcm.construct_minmover()
        else:
            minmover = BoundedMinMover(self.config.syminfo.cubicboundary, self.scorefxn)
        minimize_trial = TrialMover(minmover, dock_mcm.mc)
        minimize_trial.apply(pose)
        self.first_min_time.append(time.time() - start)
        print("first_min = time:", np.array(self.first_min_time).mean())

        # 3. docking cycle with packing and potential minimization 8 times (repack_times), accept each move with MC
        start = time.time()
        if self.config.old_behavior:
            for cycle1 in range(self.num_of_first_cycle):
                dock_mcm.apply(pose)
            # used to be a filter here in the C++ code but is now commented out.
            dock_mcm.reset_cycle_index()  # sets next_move_ to 0
            for cycle2 in range(self.num_of_second_cycle):
                dock_mcm.apply(pose)
        else:
            dock_mcm.apply(pose)
        self.dock_mcm_time.append(time.time() - start)
        print("dock_mcm = time:", np.array(self.dock_mcm_time).mean())

        # 4. Minimization, accept with MC
        # try to minimize the pose. If the pose is within bounds then accept else not
        start = time.time()
        minimize_trial.apply(pose)
        self.final_min_time.append(time.time() - start)
        print("last_min = time:", np.array(self.final_min_time).mean())

        # retrieve the lowest E pose found
        dock_mcm.mc.recover_low(pose)

class DockNRelaxProtocol:

    def __init__(self, pose, config, dofspecification: dict = default_HF_dofs):
        self.min_option_bb = config.min_option_bb
        self.min_option_sc = config.min_option_sc
        self.tf = self.create_taskfactory(pose)
        self.mf = self.create_interface_mf()
        self.cartesian = config.cartesian
        raise NotImplementedError
        self.scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
        if self.cartesian:
            self.scorefxn_cart = ScoreFunctionFactory.create_score_function("ref2015_cart")
        self.dofspecification = dofspecification

    def apply(self, pose):
        """Applies docking followed by relax and accepts according to MC."""
        start = time.time()
        docknrelax = self.create_protocol(pose)
        docknrelax.apply(pose)
        docknrelax.mc().recover_low(pose)
        print("Docking time:", time.time() - start, "seconds")

    def create_fastrelax(self):
        """Creates the fastrelax protocol."""
        fastrelax = FastRelax(1)
        fastrelax.set_scorefxn(self.scorefxn)
        fastrelax.set_task_factory(self.tf)
        fastrelax.set_movemap_factory(self.mf)
        return fastrelax

    # COPIED EXACTLY FROM ABOVE
    def create_rigidbodymover(self, pose):
        """Construct an alternative to the RigidBodyPerturbMover."""
        rb_mover = RigidBodyDofAdaptiveMover("all")
        for jump_name, jumpdof_params in self.dofspecification.items():
            for dof_name, dof_params in jumpdof_params.items():
                rb_mover.add_jump(pose, jump_name, dof_name, *self.get_extra_options(dof_params))
        return rb_mover

    # COPIED EXACTLY FROM ABOVE
    def get_extra_options(self, dof_params):
        """HACK: You cannot parse individual keyword parameters to these c++ objects, so you
        have specify them all and then change them before parsing them as below"""
        if dof_params: # check for not empty
            default = {"step_type": "gauss", "param1": 0.5, "param2": 0.0, "min_pertubation": 0.01,
                               "limit_movement": False, "max": 0, "min": 0}
            default.update(dof_params)
            return default.values()
        else:
            return []

    def create_protocol(self, pose):
        """Creates the entire dock and relax protocol."""
        dock_n_relax = SequenceMover()
        dock_n_relax.add_mover(self.create_rigidbodymover(pose))
        dock_n_relax.add_mover(self.create_fastrelax())
        return TrialMover(dock_n_relax, MonteCarlo(pose, self.scorefxn, 0.8))

    # COPIED EXACTLY FROM ABOVE
    def create_interface_mf(self):
        # construct mf and set all jumps to true
        mf = MoveMapFactory()
        mf.all_jumps(True)
        # create residueselector that just selects the interface
        residueselector = TaskSelector()
        residueselector.set_select_designable(True)
        residueselector.set_select_packable(True)
        residueselector.set_select_fixed(False)
        residueselector.set_task_factory(self.tf)
        if self.min_option_bb:
            mf.add_bb_action(move_map_action.mm_enable, residueselector)
        if self.min_option_sc:
            mf.add_chi_action(move_map_action.mm_enable, residueselector)
        # just make sure no bondangles/lengths are minimized
        mf.all_bondangles(False)
        mf.all_bondlengths(False)
        return mf


    # COPIED EXACTLY FROM ABOVE
    def create_taskfactory(self, pose, cb_dist_cutoff: int = 10.0, nearby_atom_cutoff: int = 5.5,
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
        self.cubicboundary = CubicBoundary(self.input_symdef, pose_at_initial_position=pose, dof_spec=self.dof_spec, sd=self.bound_penalty)
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