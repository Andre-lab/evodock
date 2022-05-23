#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Symmetry classes: SymDockMCMProtocol and SymDockingSlideIntoContact
@Author: Mads Jeppesen
@Date: 6/30/21
"""

from pyrosetta.rosetta.protocols.moves import MonteCarlo, TrialMover, SequenceMover, JumpOutMover, CycleMover
from pyrosetta.rosetta.protocols.rigid import RigidBodyDofSeqPerturbMover
from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMover, MinMover, PackRotamersMover
# from pyrosetta.rosetta.protocols.docking import SidechainMinmover, DockMinMover
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.task_operations import RestrictToInterfaceVectorOperation
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, InitializeFromCommandline, IncludeCurrent,\
    NoRepackDisulfides, ReadResfile, AppendRotamerSet
from pyrosetta.rosetta.core.pack.rotamer_set import UnboundRotamersOperation
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.utility import vector1_unsigned_long
from pyrosetta.rosetta.protocols.symmetry import SymDockingSlideIntoContact
from pyrosetta import SwitchResidueTypeSetMover
from src.utils import get_translation, get_rotation_euler
from pyrosetta.rosetta.core.conformation.symmetry import SlideCriteriaType
from pyrosetta.rosetta.protocols.symmetry import SequentialSymmetrySlider

class SequentialSymmetrySliderWrapper:
    """Does the the as SequentialSymmetrySlider but it needs to be constructed on every single pose to actually do
    sliding. If you just use the same SequentialSymmetrySlider on different poses it will only slide on the first
    pose."""

    def __init__(self):
        self.FA_REP_SCORE = SlideCriteriaType(2)

    def apply(self, pose):
        SequentialSymmetrySlider(pose, self.FA_REP_SCORE).apply(pose)

class SymDockingSlideIntoContactWrapper2:
    """Does the same as SymDockingSlideIntoContactWrapper but changes the underlying mover. The
    SymDockingSlideIntoContact fails if too many internal moves (1000) have been applied. We dont want out protocol
    to fail like this or we are going to have a bad time. Therefore apply of SymDockingSlideIntoContact is kept the same
    but the aspect that is calling utility_exit() is removed."""

    def __init__(self, pose):
        """Initializes the mover."""
        # this constructor needs to be called or rosetta will segfault. Reason is that the score function
        # internally defined in the mover is not initialized properly.
        dofs = pose.conformation().Symmetry_Info().get_dofs()
        self.symdockingslideintocontact = SymDockingSlideIntoContact(dofs)
        # self.
        self.switch_to_centroid = SwitchResidueTypeSetMover("centroid")

    def apply(self, pose):
        """See class description. Dockslides the pose."""
        centroid_pose = pose.clone()
        self.switch_to_centroid.apply(centroid_pose)
        # todo: this can give an error if to many slide attemps have been attempted!
        self.symdockingslideintocontact.apply(centroid_pose)
        # apply the new position of centroid_pose to pose
        self.apply_conformation_from_pose(centroid_pose, pose)

    def apply_conformation_from_pose(self, pose_from, pose_to):
        """Applies the movable dofs from pose_from to pose_to"""
        dofs = pose_from.conformation().Symmetry_Info().get_dofs()
        for jump_id, symdof in dofs.items():
            pose_to.set_jump(jump_id, pose_from.jump(jump_id))

class SymDockingSlideIntoContactWrapper:
    """The original SymDockingSlideIntoContact gives a 'Fatal Error in VDW_Energy' ERROR, which DockingSlideIntoContact"
    does not. This is because it is implemented wrong in the c++ code. in DockingSlideIntoContact the scorefxn is
    overwritten to interchain_cen which might be the problem.

    This mover attempts to overwrite the SymDockingSlideIntoContact mover so that before calling its internal apply,
    it creates pose in centroid mode that then gets supplied to the SymDockingSlideIntoContact mover. This dofs are
    then extracted from the centroid pose and supplied to the original pose."""

    def __init__(self, pose):
        """Initializes the mover."""
        # this constructor needs to be called or rosetta will segfault. Reason is that the score function
        # internally defined in the mover is not initialized properly.
        dofs = pose.conformation().Symmetry_Info().get_dofs()
        self.symdockingslideintocontact = SymDockingSlideIntoContact(dofs)
        # self.
        self.switch_to_centroid = SwitchResidueTypeSetMover("centroid")

    def apply(self, pose):
        """See class description. Dockslides the pose."""
        centroid_pose = pose.clone()
        self.switch_to_centroid.apply(centroid_pose)
        # todo: this can give an error if to many slide attemps have been attempted!
        self.symdockingslideintocontact.apply(centroid_pose)
        # apply the new position of centroid_pose to pose
        self.apply_conformation_from_pose(centroid_pose, pose)

    def apply_conformation_from_pose(self, pose_from, pose_to):
        """Applies the movable dofs from pose_from to pose_to"""
        dofs = pose_from.conformation().Symmetry_Info().get_dofs()
        for jump_id, symdof in dofs.items():
            pose_to.set_jump(jump_id, pose_from.jump(jump_id))

class SymDockingSlideIntoContactWrapper:
    """The original SymDockingSlideIntoContact gives a 'Fatal Error in VDW_Energy' ERROR, which DockingSlideIntoContact"
    does not. This is because it is implemented wrong in the c++ code. in DockingSlideIntoContact the scorefxn is
    overwritten to interchain_cen which might be the problem.

    This mover attempts to overwrite the SymDockingSlideIntoContact mover so that before calling its internal apply,
    it creates pose in centroid mode that then gets supplied to the SymDockingSlideIntoContact mover. This dofs are
    then extracted from the centroid pose and supplied to the original pose."""

    def __init__(self, pose):
        """Initializes the mover."""
        # this constructor needs to be called or rosetta will segfault. Reason is that the score function
        # internally defined in the mover is not initialized properly.
        dofs = pose.conformation().Symmetry_Info().get_dofs()
        self.symdockingslideintocontact = SymDockingSlideIntoContact(dofs)
        # self.
        self.switch_to_centroid = SwitchResidueTypeSetMover("centroid")

    def apply(self, pose):
        """See class description. Dockslides the pose."""
        centroid_pose = pose.clone()
        self.switch_to_centroid.apply(centroid_pose)
        # todo: this can give an error if to many slide attemps have been attempted!
        self.symdockingslideintocontact.apply(centroid_pose)
        # apply the new position of centroid_pose to pose
        self.apply_conformation_from_pose(centroid_pose, pose)

    def apply_conformation_from_pose(self, pose_from, pose_to):
        """Applies the movable dofs from pose_from to pose_to"""
        dofs = pose_from.conformation().Symmetry_Info().get_dofs()
        for jump_id, symdof in dofs.items():
            pose_to.set_jump(jump_id, pose_from.jump(jump_id))

class SymDockMCMProtocol:
    """To my knowledge there's no equivalent symmetrical version of DockMCMProtocol in Rosetta. This class is an
    attempt at reproducing the DockMCMProtocol in a symmetrical form and as it is used in EvoDOCK. The following
    problems are:

    1. The use of the RigidBodyPerturbMover.
       It does not consider symmetry and so will mess it up if used.
    2. The use of DockingNoRepack.
       It seems like the  mover also doesnt work for symmetry. It defines only a single jump.

    The protocol seems to do the following:
    5+45 outer protocol of:
        1. 8 x inner protocol of:
            1. Apply RigidBodyPerturbMover and RotamerTrialsMover and if E < 15.0 apply a MinMover, and accept by MCM.
            2. Same as 1 but also apply a PackRotamersMover and accept by MCM.
    """
    def __init__(self, pose, num_of_first_cycle: int = 5, num_of_second_cycle: int = 45):
        """Initializes the mover.

        :param pose: Pose to dock.
        :param num_of_first_cycle: Number of iterations for the first cycle
        :param num_of_second_cycle: Number of iterations for the second cycle
        """
        self.num_of_first_cycle = num_of_first_cycle
        self.num_of_second_cycle = num_of_second_cycle
        self.scorefxn = ScoreFunctionFactory.create_score_function("ref2015")
        self.mc = MonteCarlo(self.scorefxn, 0.8)
        self.taskkfactory = self.create_taskfactory(pose)
        self.dock_mcm_cycle = self.setup_protocol(pose)
        self.dockminmover = self.create_dockminmover()

    def apply(self, pose):
        """Applies the mover."""
        if self.mc.last_accepted_pose().empty():
            self.mc.reset(pose)
        # initial minimization
        self.dockminmover.apply(pose)
        # docking
        for _ in range(self.num_of_first_cycle):
            self.dock_mcm_cycle.apply(pose)
        # used to be a filter here in the c++ code but is now commented out.
        self.dock_mcm_cycle.reset_cycle_index() # sets next_move_ to 0
        for _ in range(self.num_of_second_cycle):
            self.dock_mcm_cycle.apply(pose)
        # minization
        self.dockminmover.mover().tolerance(0.01) #
        self.dockminmover.apply(pose)
        # recover the lowest energy pose
        self.mc.recover_low(pose)

    def create_dockminmover(self):
        """Creates a TrialMover that attemps to mimick DockMinMover."""
        min_mover = MinMover()  # movemap is initialzied
        min_mover.score_function(self.scorefxn)
        return TrialMover(min_mover, self.mc)

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
        # -- create RestrictToInterface -- #
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
        # return it
        return taskfactory

    def setup_protocol(self, pose, repack_period = 8, minimization_threshold: int = 15.0, rot_mag: int = 1.0,
                       trans_mag: int = 3.0) -> SequenceMover:
        """Sets up the protocol."""
        # todo assert pose is symmetric
        # -- Construct the meat of the dock_mcm_cycle -- #
        # construct first part
        rb_pack_min = SequenceMover()
        dofs = pose.conformation().Symmetry_Info().get_dofs()
        rb_pack_min.add_mover(RigidBodyDofSeqPerturbMover(dofs, rot_mag, trans_mag))
        rb_pack_min.add_mover(RotamerTrialsMover(self.scorefxn, self.taskkfactory))
        min_mover = MinMover() # movemap is initialzied
        min_mover.score_function(self.scorefxn) # i believe it is ref2015 by default
        rb_mover_min = JumpOutMover(rb_pack_min, min_mover, self.scorefxn, minimization_threshold)
        rb_mover_min_trial = TrialMover(rb_mover_min, self.mc)
        # construct second part
        repack_step = SequenceMover()
        repack_step.add_mover(rb_mover_min_trial)
        packer = PackRotamersMover(self.scorefxn)
        packer.task_factory(self.taskkfactory)
        repack_step.add_mover(TrialMover(packer, self.mc))
        # -- add the first and second part to the dock_mcm_cycle a repack_period amound of time -- #
        dock_mcm_cycle = CycleMover()
        for i in range(repack_period):
            dock_mcm_cycle.add_mover(rb_mover_min_trial)
            dock_mcm_cycle.add_mover(repack_step)
        return dock_mcm_cycle
