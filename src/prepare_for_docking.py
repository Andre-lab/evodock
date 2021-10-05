#!/usr/bin/env python
# coding: utf-8


from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.core.chemical import aa_vrt
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.core.pose import (addVirtualResAsRoot, chain_end_res,
                                         get_center_of_mass,
                                         remove_virtual_residues,
                                         virtual_type_for_pose)
from pyrosetta.rosetta.protocols.docking import setup_foldtree


def addVirtualResAtPoint(xyz, pose, r_start, r_end):
    d_min = 1000
    this_d = 1
    for i in range(r_start, r_end):
        rsd = pose.residue(i)
        if rsd.aa() == aa_vrt:
            continue
        if not rsd.is_protein():
            continue

        atom = rsd.atom("CA")
        this_d = (atom.xyz() - xyz).length()
        if this_d < d_min:
            d_min = this_d
            i_min = i

    rsd_type = virtual_type_for_pose(pose)
    new_res = ResidueFactory().create_residue(rsd_type)

    # move to <xyz>
    for j in range(1, new_res.natoms() + 1):
        new_res.atom(j).xyz(new_res.atom(j).xyz() + xyz)

    i_min = 2
    pose.append_residue_by_jump(new_res, i_min)

    # update PDBinfo
    pose.pdb_info().chain(pose.size(), "w")
    pose.pdb_info().number(pose.size(), 1)
    pose.pdb_info().obsolete(False)
    return pose


class PrepareForDocking:
    def __init__(self, pose):
        pose1 = Pose(pose, 1, chain_end_res(pose, 1))
        pose2 = Pose(pose, chain_end_res(pose, 1) + 1, pose.total_residue())
        pose1.pdb_info().name("chainA")
        pose2.pdb_info().name("chainB")
        self.com1 = get_center_of_mass(pose1)
        self.com2 = get_center_of_mass(pose2)

    def apply(self, pose):
        remove_virtual_residues(pose)
        addVirtualResAsRoot(self.com1, pose)
        pose = addVirtualResAtPoint(
            self.com2, pose, chain_end_res(pose, 1) + 1, pose.total_residue()
        )
        setup_foldtree(pose, "A_B", Vector1([1]))
        return pose
