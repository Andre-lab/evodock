#!/usr/bin/env python
# coding: utf-8


from src.distance_axes import calculate_local_coordinates
from src.utils import convert_range, get_position_info
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from src.utils import get_rotation_euler, get_translation

class GlobalGenotypeConverter:
    def __init__(self, native_pose, max_trans=70, syminfo: dict = None):
        if is_symmetric(native_pose):
            self.bounds = self.define_symmetric_bounds(native_pose, syminfo)
        else:
            self.max_rot = 180
            mtrans = max_trans
            self.max_trans = [mtrans for t in range(3)]
            self.min_trans = [mtrans * -1 for t in range(3)]
            self.bounds = self.define_bounds()

    def define_symmetric_bounds(self, native_pose, syminfo):
        """Define symmetrical bounds."""
        bounds = []
        for jump, dofs, parsed_bounds in zip(syminfo.get("jumps_int"), syminfo.get("dofs_int"), syminfo.get("bounds")):
            flexible_jump = native_pose.jump(jump)
            rot = get_rotation_euler(flexible_jump)
            trans = get_translation(flexible_jump)
            for dof, bound in zip(dofs, parsed_bounds):
                if dof < 4:
                    native_val = trans[dof - 1]
                else:
                    native_val = rot[dof - 4]
                bounds.append((- float(bound) + native_val, + float(bound) + native_val))
        return bounds

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
            bounds.append(
                (
                    init_pos_trans[i] - self.min_trans[i],
                    init_pos_trans[i] - self.max_trans[i],
                )
            )
        return bounds


def main():
    print("genotype converter")


if __name__ == "__main__":
    main()
