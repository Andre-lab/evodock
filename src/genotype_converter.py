#!/usr/bin/env python
# coding: utf-8


from src.distance_axes import calculate_local_coordinates
from src.utils import convert_range, get_position_info


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
        self.max_rot = [3, 3, 3]
        self.max_trans = [1, 1, 1]
        self.pose = native_pose
        self.bounds = self.define_bounds()

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
