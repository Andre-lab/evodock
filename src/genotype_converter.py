#!/usr/bin/env python
# coding: utf-8

import numpy as np

from src.distance_axes import (calculate_local_coordinates,
                               calculate_max_coordiantes)
# from principal_axes import calculate_max_coordiantes
from src.utils import convert_range, get_position_info


def get_translation_max(dock_pose):
    jump_num = 1
    flexible_jump = dock_pose.jump(jump_num)
    translation = np.asarray(flexible_jump.get_translation())
    max_trans = max(translation) + 1
    return [max_trans, max_trans, max_trans]


class GlobalGenotypeConverter:
    def __init__(self, native_pose):
        self.max_rot = 180
        # self.max_trans = get_translation_max(native_pose)
        max_trans_coordinates, min_trans_coordinates = calculate_max_coordiantes(
            native_pose
        )
        self.max_trans = [t + 10 for t in max_trans_coordinates]
        self.min_trans = [t - 10 for t in min_trans_coordinates]
        # self.max_trans = [70, 70, 70]
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
        self.max_rot = 16
        self.pose = native_pose
        self.max_trans, self.min_trans = calculate_local_coordinates(native_pose)
        self.bounds = self.define_bounds()

    def define_bounds(self):
        bounds = []
        init_pos = get_position_info(self.pose)
        init_pos_rot = init_pos[:3]
        init_pos_trans = init_pos[3:]
        for i in range(3):
            bounds.append(
                (init_pos_rot[i] - self.max_rot, init_pos_rot[i] + self.max_rot)
            )
        for i in range(3):
            bounds.append((self.min_trans[i], self.max_trans[i],))
        return bounds


def main():
    print("genotype converter")


if __name__ == "__main__":
    main()
