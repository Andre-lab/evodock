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
        self.size = len(self.bounds)

    # Handles both local search and global search cases!!
    # For local docking: Use -initialize_rigid_body_dof and set the bounds to whatever you want to search
    # For global docking: Set the bounds very high with and with -initialize_rigid_body_dof, with it give the starting
    #   position to be that including the bounds you give it.
    # if you dont parse -initialize_rigid_body_dof the native_val will be 0, otherwise the value set in the symdef file.
    def define_symmetric_bounds(self, native_pose, syminfo):
        """Define symmetrical bounds."""
        bounds = []
        for jump, dofs, parsed_bounds in zip(syminfo.get("jumps_int"), syminfo.get("dofs_int"), syminfo.get("bounds")):
            flexible_jump = native_pose.jump(jump)
            rot = get_rotation_euler(flexible_jump)
            trans = get_translation(flexible_jump)
            for dof, bound in zip(dofs, parsed_bounds):
                # rotation is from maximum from -180 to +180
                # trans is from 0->inf and shouldnt go into the negatives (can mess up the symmetry)
                if dof < 4: # translational dof
                    native_val = trans[dof - 1]
                    min_value = max(0, native_val - float(bound)) # make sure the lowest value is 0
                    max_value = native_val + float(bound)
                    bounds.append((min_value, max_value))
                else: # rotational dof
                    native_val = rot[dof - 4]
                    min_value = max(-180, native_val - float(bound)) # makes sure the lowest value is -180
                    max_value = min(180,  native_val + float(bound)) # makes sure the highest value is 180
                    bounds.append((min_value,max_value))
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
