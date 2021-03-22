#!/usr/bin/env python
# coding: utf-8


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


class PosePacker:
    def __init__(self, pose, scorefxn):
        self.pose = pose
        self.scorefxn = scorefxn
        # SSRB: A repack was necessary as high fa_rep was being observed after side chain recovery
        local_tf = TaskFactory()
        local_tf.push_back(IncludeCurrent())
        local_tf.push_back(RestrictToRepacking())
        local_tf.push_back(NoRepackDisulfides())
        conformer_full_repack = PackRotamersMover(self.scorefxn)
        conformer_full_repack.task_factory(local_tf)
        self.packer = conformer_full_repack

    def run(self):
        # pack_pdb = self.pack(self.input_pose)
        # pack_pdb.dump_pdb(self.output_name)

        return pack_pdb

    def apply(self, pose):
        self.packer.apply(pose)

    def pack(self, pose):
        pack_pose = Pose()
        pack_pose.assign(pose)

        pose_packer = standard_packer_task(pack_pose)

        pose_packer.restrict_to_repacking()
        pose_packer.or_include_current(True)

        # SSRB: A repack was necessary as high fa_rep was being observed after side chain recovery
        local_tf = TaskFactory()
        local_tf.push_back(IncludeCurrent())
        local_tf.push_back(RestrictToRepacking())
        local_tf.push_back(NoRepackDisulfides())
        conformer_full_repack = PackRotamersMover(self.scorefxn)
        conformer_full_repack.task_factory(local_tf)

        # ptask = PackerTask()
        # ptask.include_current(
        # packmover = PackRotamersMover(self.scorefxn, pose_packer)
        # packmover.apply(pack_pose)
        # return pack_pose
        return conformer_full_repack


class CustomRotamer:
    def __init__(self, pose, scorefxn_):
        # setup the movemap
        # tf_ = TaskFactory()
        # tf_.push_back(IncludeCurrent())
        # tf_.push_back(RestrictToRepacking())
        # tf_.push_back(NoRepackDisulfides())
        mcm_docking = DockMCMProtocol()
        mcm_docking.set_native_pose(pose)
        mcm_docking.set_scorefxn(scorefxn_)
        mcm_docking.set_rt_min(False)
        mcm_docking.set_sc_min(False)
        mcm_docking.apply(pose)
        # rotamer trial
        self.packer = RotamerTrialsMover(scorefxn_, mcm_docking.task_factory())

    def apply(self, pose):
        self.packer.apply(pose)


class CustomPacker:
    def __init__(self, pose, scorefxn_):
        # pose packer
        # tf_ = TaskFactory()
        # pose_packer = standard_packer_task(pose)
        # pose_packer.restrict_to_repacking()
        # pose_packer.or_include_current(True)
        # self.packmover = PackRotamersMover(scorefxn_, pose_packer)
        # mcm_docking = DockMCMProtocol()
        # mcm_docking.set_native_pose(pose)
        # mcm_docking.set_scorefxn(scorefxn_)
        # mcm_docking.set_rt_min(False)
        # mcm_docking.set_sc_min(False)
        # mcm_docking.apply(pose)
        # print("FINISH !! --------------------------------")
        # self.packmover = PackRotamersMover(scorefxn_)
        # self.packmover.task_factory(mcm_docking.task_factory())

        # SSRB: A repack was necessary as high fa_rep was being observed after side chain recovery
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
        conformer_full_repack = PackRotamersMover(scorefxn_)
        conformer_full_repack.task_factory(local_tf)
        self.packmover = conformer_full_repack

    def apply(self, pose):
        self.packmover.apply(pose)


class CustomPackRotamer:
    def __init__(self, pose, scorefxn_):
        # setup the movemap
        movemap_ = MoveMap()
        movemap_.set_chi(False)
        movemap_.set_bb(False)
        movemap_.set_jump(1, True)

        min_tolerance_ = 0.01
        trans_magnitude_ = 0.1
        rot_magnitude_ = 5.0
        min_type_ = "lbfgs_armijo_nonmonotone"
        nb_list_ = True

        minimization_threshold = 15.0
        mc_ = MonteCarlo(scorefxn_, 0.8)
        tf_ = TaskFactory()

        rb_mover = RigidBodyPerturbMover(
            pose, movemap_, rot_magnitude_, trans_magnitude_, partner_downstream, True
        )

        # rotamer trial
        self.packer = RotamerTrialsMover(scorefxn_, tf_)

        # pose packer
        pose_packer = standard_packer_task(pose)
        pose_packer.restrict_to_repacking()
        pose_packer.or_include_current(True)
        self.packmover = PackRotamersMover(scorefxn_, pose_packer)

        rb_pack_min = SequenceMover()
        rb_pack_min.add_mover(rb_mover)
        rb_pack_min.add_mover(self.packer)

        min_mover = MinMover(movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_)
        rb_mover_min = JumpOutMover(
            rb_pack_min, min_mover, scorefxn_, minimization_threshold
        )

        self.rb_mover_min_trail = TrialMover(rb_mover_min, mc_)

        self.repack_step = SequenceMover()
        self.repack_step.add_mover(self.rb_mover_min_trail)
        pack_rotamers = PackRotamersMover(scorefxn_)
        self.pack_interface_and_move_loops_trial = TrialMover(pack_rotamers, mc_)
        self.repack_step.add_mover(self.pack_interface_and_move_loops_trial)

    def apply(self, pose):
        self.repack_step.apply(pose)


class RigidBodyMover:
    def __init__(self, pose, scorefxn_):
        # setup the movemap
        movemap_ = MoveMap()
        movemap_.set_chi(False)
        movemap_.set_bb(False)
        movemap_.set_jump(1, True)

        min_tolerance_ = 0.01
        trans_magnitude_ = 0.1
        rot_magnitude_ = 5.0
        min_type_ = "lbfgs_armijo_nonmonotone"
        nb_list_ = True

        minimization_threshold = 15.0
        mc_ = MonteCarlo(scorefxn_, 0.8)
        tf_ = TaskFactory()

        rb_mover = RigidBodyPerturbMover(
            pose, movemap_, rot_magnitude_, trans_magnitude_, partner_downstream, True
        )

        # rotamer trial
        self.packer = RotamerTrialsMover(scorefxn_, tf_)

        # pose packer
        pose_packer = standard_packer_task(pose)
        pose_packer.restrict_to_repacking()
        pose_packer.or_include_current(True)
        self.packmover = PackRotamersMover(scorefxn_, pose_packer)

        # packer = PackRotamersMover(scorefxn_)
        # packer = PosePacker(pose, scorefxn_)

        rb_pack_min = SequenceMover()
        rb_pack_min.add_mover(rb_mover)
        rb_pack_min.add_mover(self.packer)

        min_mover = MinMover(movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_)
        rb_mover_min = JumpOutMover(
            rb_pack_min, min_mover, scorefxn_, minimization_threshold
        )

        self.rb_mover_min_trail = TrialMover(rb_mover_min, mc_)

        self.repack_step = SequenceMover()
        self.repack_step.add_mover(self.rb_mover_min_trail)
        pack_rotamers = PackRotamersMover(scorefxn_)
        self.pack_interface_and_move_loops_trial = TrialMover(pack_rotamers, mc_)
        self.repack_step.add_mover(self.pack_interface_and_move_loops_trial)

    def apply(self, pose):
        self.packmover.apply(pose)
