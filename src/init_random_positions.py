#!/usr/bin/env python
# coding: utf-8


from pyrosetta import Pose

from src.perturbation_step import PerturbationStepMover, conversion_movers
from src.reader import StructureReader

pdb_archive = {
    "1oph": "./easy_dock/1oph_AB.prepack.pdb",
    "1ml0": "./easy_dock/1ml0_AD.prepack.pdb",
    "1ktz": "./easy_dock/1ktz_AB.prepack.pdb",
    "1qa9": "./easy_dock/1qa9_AB.prepack.pdb",
    # "2hrk": "./easy_dock/2hrk_AB_0001.prepack.pdb",
    "2hrk": "./easy_dock/2hrk_AB.prepack.pdb",
    "1kxp": "./easy_dock/1kxp_AD.prepack.pdb",
    "2hle": "./easy_dock/2hle_AB.prepack.pdb",
    "1b6c": "./easy_dock/1b6c_AB.prepack.pdb",
    "1ppe": "./easy_dock/1ppe_IE.prepack.pdb",
}


def start_input_poses(pose_input):
    pose = StructureReader().get_from_file(pose_input)
    prot = [p for p in pdb_archive.keys() if p in pose_input][0]
    native_pose = StructureReader().get_from_file(pdb_archive[prot])

    # 3. create centroid <--> fullatom conversion Movers
    to_centroid, to_fullatom, recover_sidechains = conversion_movers(pose)

    # perturb = PerturbationStepMover().initialize(
    #     pose=pose,
    #     dock_jump=1,
    #     translation=12,
    #     rotation=35,
    #     to_fullatom=to_fullatom,
    #     recover_sidechains=recover_sidechains,
    # )

    # perturb.apply(pose)
    # to_centroid.apply(pose)
    return native_pose, pose
