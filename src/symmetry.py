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

from pyrosetta.rosetta.protocols.moves import MonteCarlo, TrialMover, SequenceMover, JumpOutMover, CycleMover
from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMover, MinMover, PackRotamersMover
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.task_operations import RestrictToInterfaceVectorOperation
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, InitializeFromCommandline, IncludeCurrent,\
    NoRepackDisulfides, ReadResfile, AppendRotamerSet
from pyrosetta.rosetta.core.pack.rotamer_set import UnboundRotamersOperation
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.utility import vector1_unsigned_long
from pyrosetta.rosetta.protocols.rigid import RigidBodyDofAdaptiveMover
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_jump_num, is_symmetric, jump_num_sym_dof
from pyrosetta.rosetta.core.kinematics import MoveMap

default_dofs = {"JUMPHFfold1": {
                    "z": {"param1": 0.1},
                    "z_angle": {"param1": 1} # For the future{"limit_movement": True},
                },
                "JUMPHFfold111": {
                    "x": {"param1": 0.1}
                },
                "JUMPHFfold1111": {
                    "x_angle": {"param1": 3},
                    "y_angle": {"param1": 3},
                    "z_angle": {"param1": 3}}}

class SymDockMCMCycle:
    """The symmetric version of DockMCMCycler class in C++ with extra functionality. The extra functionality comes from the posibility
    of specifying a dofspecification that fine tunes seperate symmetric dofs. See the comments in setup_protocol for a
    description of the algorithm."""

    def __init__(self, pose, repack_period = 8, dofspecification: dict = default_dofs):
        """Initialization"""
        self.set_default(pose, repack_period, dofspecification)
        self.dock_mcm_cycle = None

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

    def set_default(self, pose, repack_period, dof_specification, bb_min_res:list=None, sc_min_res:list=False):
        """Sets up the the default parameteres used when creating the dock_mcm_cycle object in the 'setup_protocol' function."""
        self.dofspecification = dof_specification
        self.sc_min = False
        self.rt_min = False
        self.scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
        self.movemap = MoveMap()
        self.movemap.set_chi(False)
        self.movemap.set_bb(False)
        for jump_name, jumpdof_params in self.dofspecification.items():
            jumpid = sym_dof_jump_num(pose, jump_name)
            self.movemap.set_jump(jumpid, True)
        # TODO: check if these are defaulted to all residues?
        if bb_min_res or sc_min_res:
            raise NotImplementedError()
        self.min_tolerance = 0.01
        self.min_type = "lbfgs_armijo_nonmonotone"
        self.nb_list = True

        # DockMCM mover only initializes an empty taskfatory so one have to set it manually here
        self.tf = self.create_taskfactory(pose)
        self.mc = MonteCarlo(self.scorefxn, 0.8)
        self.repack_period = repack_period

    def setup_protocol(self, pose):
        """Sets up the dock_mcm_cycle object."""
        # Contruction of rb_mover_min_trial (The first move in the sequence)
        # 1: (1a. rb dock move -> 1b. rotamer trial) -> 2: if E < 15.0 apply minimization.
        # Accept the entire move 1+2 with the MC criterion
        self.rb_mover = self.construct_rigidbodymover(pose)
        rottrial = RotamerTrialsMover(self.scorefxn, self.tf)
        rb_pack_min = SequenceMover()
        rb_pack_min.add_mover(self.rb_mover)
        rb_pack_min.add_mover(rottrial)
        minimization_threshold = 15.0
        min_mover = MinMover(self.movemap, self.scorefxn, self.min_type, self.min_tolerance, self.nb_list)
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
        if self.rt_min or self.sc_min:
            raise NotImplementedError()

        # Putting the 2 moves (rb_mover_min_trial and repack_step) into a CycleMover that is repeat self.repack_period times
        dock_mcm_cycle = CycleMover()
        for i in range(self.repack_period):
            dock_mcm_cycle.add_mover(rb_mover_min_trial)
            dock_mcm_cycle.add_mover(repack_step)

        self.dock_mcm_cycle = dock_mcm_cycle

    def construct_rigidbodymover(self, pose):
        """Construct an alternative to the RigidBodyPerturbMover."""
        rb_mover = RigidBodyDofAdaptiveMover("all")
        for jump_name, jumpdof_params in self.dofspecification.items():
            for dof_name, dof_params in jumpdof_params.items():
                rb_mover.add_jump(pose, jump_name, dof_name, *self.get_extra_options(dof_params))
        return rb_mover

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

class SymDockMCMProtocol:
    """The symmetric version of DockMCMCycler class in C++ with extra functionality.

    The Algorithm is as follows:
    1. Fixed bb packing (PackRotamersMover) - accept by MC.
    2.
    """

    # The protocol seems to do the following:
    # <num_of_first_cycle>+<num_of_second_cycle> outer protocol of:
    #     1. 8 x inner protocol of:
    #         1. Apply RigidBodyPerturbMover and RotamerTrialsMover and if E < 15.0 apply a MinMover, and accept by MCM.
    #         2. Same as 1 but also apply a PackRotamersMover and accept by MCM.

    def __init__(self, num_of_first_cycle: int = 4, num_of_second_cycle: int = 45):
        """Initializes the mover.

        :param pose: Pose to dock.
        :param num_of_first_cycle: Number of iterations for the first cycle
        :param num_of_second_cycle: Number of iterations for the second cycle
        """
        self.num_of_first_cycle = num_of_first_cycle
        self.num_of_second_cycle = num_of_second_cycle
        self.scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
        # these are form DockingHighRes - they are currently not used in the protocol
        self.sc_min = False
        self.rt_min = False

    def apply(self, pose):
        """Applies the mover."""
        assert is_symmetric(pose), "Pose is not symmetric!"

        # 1. Packing of pose, accept with MC
        dock_mcm = SymDockMCMCycle(pose)
        # we need to call at least reset in order to initialize the pose in the MC object.
        # Since we just initialized the dock_mcm objet here it wil allways be empty. This is how
        # it is in the code but it seems a bit redundant because setup_protocol will then be called twice which is unnecessary
        if dock_mcm.mc.last_accepted_pose().empty():
            dock_mcm.init_mc(pose)
        initial_pack = PackRotamersMover()
        initial_pack.score_function(self.scorefxn)
        initial_pack.task_factory(dock_mcm.tf)
        initial_pack_trial = TrialMover(initial_pack, dock_mcm.mc)
        initial_repack_sequence = SequenceMover()
        initial_repack_sequence.add_mover(initial_pack_trial)
        if self.rt_min or self.sc_min:
            raise NotImplementedError()
        initial_repack_sequence.apply(pose)

        # 2. Minimization, accept with MC
        # Called DockMinMover in C++, but this just seems to be a minmover in a trialmover with the exact same minmover settings as dock_mcm
        minmover = MinMover(dock_mcm.movemap, dock_mcm.scorefxn, dock_mcm.min_type, dock_mcm.min_tolerance, dock_mcm.nb_list)
        minimize_trial = TrialMover(minmover, dock_mcm.mc)
        minimize_trial.apply(pose)

        # 3. docking cycle with packing and potential minimization 8 times (repack_times), accept each move with MC
        for cycle1 in range(self.num_of_first_cycle):
            dock_mcm.apply(pose)
        # used to be a filter here in the C++ code but is now commented out.
        dock_mcm.reset_cycle_index() # sets next_move_ to 0
        for cycle2 in range(self.num_of_second_cycle):
            dock_mcm.apply(pose)

        # 4. Minimization, accept with MC
        minimize_trial.apply(pose)

        # retrieve the lowest E pose found
        dock_mcm.mc.recover_low(pose)