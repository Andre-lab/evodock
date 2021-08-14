#!/usr/bin/env python
# coding: utf-8

import sys

from pyrosetta import Pose, init, pose_from_file, standard_packer_task
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent, InitializeFromCommandline, NoRepackDisulfides,
    RestrictToRepacking)
from pyrosetta.rosetta.protocols.loops import get_fa_scorefxn
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.moves import PyMOLMover


class PosePacker:
    def __init__(self, pose, filename, mode_="standard"):
        self.pose = pose
        self.output_name = filename.replace(".pdb", ".prepack.pdb")
        self.mode = mode_
        self.scorefxn = get_fa_scorefxn()

    def run(self):
        pack_pdb = self.pack()
        pack_pdb.dump_pdb(self.output_name)
        return pack_pdb

    def pack(self):
        if self.mode == "standard":
            return PoseStandardPacker(self.pose, self.output_name).pack()
        if self.mode == "custom":
            return PosePackRotamers(self.pose, self.output_name).pack()


class PoseStandardPacker(PosePacker):
    def pack(self):
        pack_pose = Pose()
        pack_pose.assign(self.pose)

        pose_packer = standard_packer_task(pack_pose)
        pose_packer.restrict_to_repacking()
        pose_packer.or_include_current(False)

        packmover = PackRotamersMover(self.scorefxn, pose_packer)
        packmover.apply(pack_pose)
        return pack_pose


class PosePackRotamers(PosePacker):
    def pack(self):
        pack_pose = Pose()
        pack_pose.assign(self.pose)

        local_tf = TaskFactory()
        local_tf.push_back(InitializeFromCommandline())
        # local_tf.push_back(IncludeCurrent())
        local_tf.push_back(RestrictToRepacking())
        local_tf.push_back(NoRepackDisulfides())
        conformer_full_repack = PackRotamersMover(self.scorefxn)
        conformer_full_repack.task_factory(local_tf)
        conformer_full_repack.apply(pack_pose)
        return pack_pose


def run_repacking(filename, opts, pymover=None):
    init(extra_options=opts)
    pose = Pose()
    pose_from_file(pose, filename)
    pose.pdb_info().name("INIT_STATE")
    if pymover is not None:
        pymover.apply(pose)
    scorefxn = get_fa_scorefxn()
    print("INIT ENERGY ", scorefxn.score(pose))
    packer = PosePacker(pose, filename, "custom")
    packed_pdb = packer.run()
    final_energy = scorefxn.score(packed_pdb)
    print("FINAL ENERGY ", final_energy)
    packed_pdb.pdb_info().name("PACKED")

    if pymover is not None:
        pymover.apply(packed_pdb)
    return final_energy, packed_pdb


def main():
    pymover = PyMOLMover(address="10.8.0.6", port=65000, max_packet_size=1400)
    filename = sys.argv[-1]
    opts = " ".join(
        [
            "-mute all",
            # "-unmute core.pack.pack_rotamers core.pack.dunbrack",
            "-partners A_B",
            "-ex1",
            "-ex2aro",
            "-extrachi_cutoff 1",
            # "-unboundrot {}".format(filename),
        ]
    )
    run_repacking(filename, opts, pymover)


if __name__ == "__main__":
    main()
