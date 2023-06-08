#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Aligns 2 structures
@Author: Mads Jeppesen
@Date: 5/19/22
"""
from pyrosetta import init, PyMOLMover
from pathlib import Path
from pyrosetta import init, pose_from_file
from pyrosetta.rosetta.core.conformation.symmetry import residue_center_of_mass
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.core.scoring import CA_rmsd
import argparse

def main(fix, rot, outname, outrmsd=None):
    # make starting pose
    fix_pose = pose_from_file(fix)
    # r = residue_center_of_mass(fix_pose.conformation(), 1, fix_pose.size())
    rot_pose = pose_from_file(rot)
    calpha_superimpose_pose(rot_pose, fix_pose)
    if outrmsd:
        with open(outrmsd, "w") as f:
            f.write(CA_rmsd(fix_pose, rot_pose))
    rot_pose.dump_pdb(outname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="align", description="aligns a rotating pdb file onto a fixed pdb file")
    parser.add_argument('--fix', help="pdb file to align to", required=True)
    parser.add_argument('--rot', help="pdb file to align", required=True)
    parser.add_argument('--outname', help="out name of the file", required=True)
    parser.add_argument('--outrmsd', help="out name of a file containing the rmsd")
    args = parser.parse_args()
    init()
    main(args.fix, args.rot, args.outname, args.outrmsd)