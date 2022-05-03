"""New Document"""
#!/usr/bin/env python
# coding: utf-8

import os
import glob
import sys
import argparse


from sys import argv, stderr, stdout
from pyrosetta import init, pose_from_file, Pose
from pyrosetta.rosetta.core.pose import chain_end_res
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.toolbox.mutants import mutate_residue


def get_pose_from_file(input_file):
    new_pose = Pose()
    pose_from_file(new_pose, input_file)
    return new_pose


def join_two_chunks(C1_pose, C2_pose, new_chain=True):
    new_pose = Pose()
    new_pose.assign(C1_pose)
    append_pose_to_pose(new_pose, C2_pose, new_chain)
    return new_pose


def renum_pdb(pdbname):
    lines = open(pdbname, "r").readlines()
    oldresnum = "   "
    count = 0
    outid = stdout

    for line in lines:
        line_edit = line
        if line[0:3] == "TER":
            continue

        if line_edit[0:4] == "ATOM" or line_edit[0:6] == "HETATM":
            if not (line[16] == " " or line[16] == "A"):
                continue

            resnum = line_edit[23:26]
            if not resnum == oldresnum:
                count = count + 1
            oldresnum = resnum

            newnum = "%3d" % count
            line_edit = line_edit[0:23] + newnum + line_edit[26:]

            outid.write(line_edit)

    outid.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed", help="path to fixed pdb")
    parser.add_argument("--rotating", help="path to rotating pdb")
    args = parser.parse_args()


if __name__ == "__main__":
    main()
