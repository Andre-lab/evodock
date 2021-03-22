#!usr/bin/env python

import logging
from math import sqrt

import numpy as np
import pyrosetta.rosetta as rosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import \
    ScoreFunction as RosettaScoreFunction
from pyrosetta.rosetta.core.scoring.shape import ZernikeDescriptorCalculator
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from pyrosetta.rosetta.protocols.toolbox.rigid_body import \
    create_euler_rotation
from scipy.spatial.transform import Rotation as R

from src.position_utils import build_axis, get_translation, to_rosetta
from src.utils import IP_ADDRESS, convert_gamma, convert_translation


def l2_norm(zd_reference, zd_current):
    L2_norm = 0
    for i, v in enumerate(zd_current):
        L2_norm += (zd_current[i] - zd_reference[i]) * (zd_current[i] - zd_reference[i])
    sqrt_L2_norm = sqrt(L2_norm)
    return sqrt_L2_norm


class FitnessFunction:
    def __init__(self, reference, native_pose, input_pose, trans_max_magnitude):
        self.logger = logging.getLogger("evodock.scfxn")
        self.native_pose = native_pose
        self.logger.setLevel(logging.INFO)
        self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.zd_calculator = ZernikeDescriptorCalculator()
        self.reference_ZD = self.zd_calculator.invariants_from_file(reference)

        self.scfxn_rosetta = RosettaScoreFunction()
        self.scfxn_rosetta.set_weight(
            rosetta.core.scoring.ScoreType.zernike_descriptor, 100
        )
        # self.reference_ZD = self.zd_calculator.invariants_from_pose(dock_pose)
        self.dock_pose = Pose()
        self.rot_max_magnitude = 180
        self.trans_max_magnitude = trans_max_magnitude

        self.pymover.apply(input_pose)
        self.dock_pose.assign(input_pose)
        self.dock_pose.pdb_info().name("INIT_STATE")
        self.pymover.apply(self.dock_pose)
        self.jump_num = 1
        self.ax1, self.ax2, self.ax3 = build_axis()
        mros_temp = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.mros, _ = to_rosetta(mros_temp, [0, 0, 0])

    def get_ZD_from_file(self, filename):
        invariants = self.zd_calculator.invariants_from_file(filename)
        return invariants

    def get_ZD_from_pose(self, pose):
        invariants = self.zd_calculator.invariants_from_pose(pose)
        return invariants

    def starting_pdb(self):
        pose = self.apply_sixD_to_pose([0, 0, 0, 0, 0, 0])
        return pose

    def get_rmsd(self, pose):
        rmsd = CA_rmsd(self.native_pose, pose)
        return rmsd

    def render_models(self, pdb_id, SixD_vector, is_best=None):
        pose = self.apply_sixD_to_pose(SixD_vector)
        dst = self.scfxn_rosetta.score(pose)
        if np.isnan(dst):
            dst = 10000
        prot_name = "popul" if is_best is None else is_best
        pose.pdb_info().name(prot_name + "_pose_" + str(pdb_id))
        self.pymover.apply(pose)
        rmsd = self.get_rmsd(pose)
        return dst, rmsd

    def score(self, SixD_vector):
        pose = self.apply_sixD_to_pose(SixD_vector)
        try:
            dst = self.scfxn_rosetta.score(pose)
        except ValueError:
            dst = 10000
        if np.isnan(dst):
            dst = 10000
        return dst

    def goal_rotation_and_translation(self):
        flexible_jump = self.dock_pose.jump(self.jump_num)
        self.logger.info("GOAL DATA: ")
        self.logger.info("init score ", self.scfxn_rosetta.score(self.dock_pose))
        self.logger.info(flexible_jump.get_rotation())
        self.logger.info(flexible_jump.get_translation())

    def size(self):
        return 6

    def convert_genotype(self, genotype):
        gen = []
        for i, g in enumerate(genotype):
            if i < 2:
                new_value = g * self.rot_max_magnitude
                gen.append(new_value)
            else:
                if i == 2:
                    gen.append(convert_gamma(g))
                else:
                    gen.append(convert_translation(g, self.trans_max_magnitude))
        return gen

    def get_sol_string(self, sol):
        return " , ".join([str(e) for e in sol])

    def apply_sixD_to_pose(self, SixD_vector):
        ind_pose = Pose()
        ind_pose.assign(self.dock_pose)
        euler = np.asarray(SixD_vector[0:3])
        r = R.from_euler("xyz", euler, degrees=True).as_matrix()
        flexible_jump = ind_pose.jump(self.jump_num)
        rosetta_rotation, rosetta_translation = to_rosetta(r, SixD_vector[3:])
        flexible_jump.set_rotation(rosetta_rotation)
        flexible_jump.set_translation(rosetta_translation)
        ind_pose.set_jump(self.jump_num, flexible_jump)
        # now is time to score the join_pose
        return ind_pose
