#!/usr/bin/env python
# coding: utf-8

import glob
import random

from pyrosetta.rosetta.core.import_pose import poses_from_files
from pyrosetta.rosetta.utility import vector1_std_string

from src.utils import convert_range, get_pose_from_file, get_position_info

from pyrosetta.rosetta.protocols.toolbox import CA_superimpose


class RefineCluspro:
    def __init__(self, config, max_trans):
        input_pose = get_pose_from_file(config.pose_input)
        lst = glob.glob(config.init_pdbs)
        filenames = vector1_std_string()
        for f in lst:
            filenames.append(f)
        self.list_models = poses_from_files(filenames)
        self.list_models = [CA_superimpose(input_pose, p) for p in self.list_models]
        self.positions = [get_position_info(p) for p in self.list_models]
        max_boundaries = [max([abs(x[i]) for x in self.positions]) for i in range(6)]
        max_boundaries = list(map(lambda i: (i, i * -1), max_boundaries))
        self.converter = GlobalGenotypeConverter(input_pose, max_trans)

    def refine_cluspro(self, idx):
        idx = int(idx % len(self.list_models))
        return self.converter.convert_positions_to_genotype(self.positions[idx])


def generate_genotype(pose_input, trans_max_magnitude):
    max_rot = [8, 8, 8]
    max_trans = [3, 3, 3]
    def set_bounds(pose_input):
        bounds = []
        init_pos = get_position_info(pose_input)
        init_rot = init_pos[:3]
        for i in range(3):
            bounds.append((init_rot[i] - max_rot[i], init_rot[i] + max_rot[i]))
        init_trans = init_pos[3:]
        for i in range(3):
            bounds.append((init_trans[i] - max_trans[i], init_trans[i] + max_trans[i]))
        return bounds

    converter = GlobalGenotypeConverter(pose_input, trans_max_magnitude)
    local_bounds = set_bounds(pose_input)
    new_values = []
    for i in range(6):
        new_values.append(random.uniform(local_bounds[i][0], local_bounds[i][1]))

    genotype = converter.convert_positions_to_genotype(new_values)
    return genotype


class GlobalGenotypeConverter:
    def __init__(self, native_pose, max_trans=70):
        self.max_rot = 180
        mtrans = max_trans
        self.max_trans = [mtrans for t in range(3)]
        self.min_trans = [mtrans * -1 for t in range(3)]
        self.bounds = self.define_bounds()

    def define_bounds(self):
        bounds = []
        # bound for rotation
        for i in range(3):
            bounds.append((self.max_rot * -1, self.max_rot))
        # bound for translation
        for i in range(3):
            bounds.append((self.min_trans[i], self.max_trans[i]))
        return bounds

    def convert_genotype(self, genotype):
        gen = []
        for i, g in enumerate(genotype):
            new_value = convert_range(
                g, (-1, 1), (self.bounds[i][0], self.bounds[i][1])
            )
            gen.append(new_value)

        return gen

    def convert_positions_to_genotype(self, positions):
        gen = []
        for i, val in enumerate(positions):
            new_value = convert_range(
                val, (self.bounds[i][0], self.bounds[i][1]), (-1, 1)
            )
            gen.append(new_value)

        return gen


class LocalGenotypeConverter(GlobalGenotypeConverter):
    def __init__(self, native_pose):
        self.max_rot = [0.01, 0.01, 0.01]
        self.max_trans = [0.01, 0.01, 0.01]
        self.pose = native_pose
        self.bounds = self.define_bounds()
        print("Local Genotype is legacy code")
        exit()

    def define_bounds(self):
        bounds = []
        init_pos = get_position_info(self.pose)
        init_rot = init_pos[:3]
        for i in range(3):
            bounds.append(
                (init_rot[i] - self.max_rot[i], init_rot[i] + self.max_rot[i])
            )
        init_trans = init_pos[3:]
        for i in range(3):
            bounds.append(
                (init_trans[i] - self.max_trans[i], init_trans[i] + self.max_trans[i])
            )
        return bounds


def main():
    print("genotype converter")


if __name__ == "__main__":
    main()
