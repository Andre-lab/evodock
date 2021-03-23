#!usr/bin/env python

import logging

import numpy as np
from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreFunctionFactory
from pyrosetta.rosetta.protocols.docking import (calc_interaction_energy,
                                                 calc_Irmsd)
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from scipy.spatial.transform import Rotation as R

from src.genotype_converter import GlobalGenotypeConverter
from src.position_utils import build_axis, to_rosetta
from src.utils import IP_ADDRESS


class FAFitnessFunction:
    def __init__(self, native_pose, trans_max_magnitude):
        self.logger = logging.getLogger("evodock.scfxn")
        self.native_pose = native_pose
        self.input_pose = Pose()
        self.input_pose.assign(self.native_pose)
        self.logger.setLevel(logging.INFO)
        self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.scfxn_rosetta = ScoreFunctionFactory.create_score_function("ref2015")
        self.dock_pose = Pose()
        self.converter = GlobalGenotypeConverter(self.input_pose)
        self.pymover.apply(self.input_pose)
        self.trans_max_magnitude = trans_max_magnitude
        self.dock_pose.assign(self.input_pose)
        self.dock_pose.pdb_info().name("INIT_STATE")
        self.pymover.apply(self.dock_pose)
        self.jump_num = 1
        self.ax1, self.ax2, self.ax3 = build_axis()
        mros_temp = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.mros, _ = to_rosetta(mros_temp, [0, 0, 0])

    def get_rmsd(self, pose):
        rmsd = CA_rmsd(self.native_pose, pose)
        return rmsd

    def render_individual(self, pdb_id, individual, is_best=None, interface=False):
        pose = individual.pose
        dst = self.scfxn_rosetta.score(pose)
        interface = calc_interaction_energy(pose, self.scfxn_rosetta, Vector1([1]))
        irms = calc_Irmsd(self.native_pose, pose, self.scfxn_rosetta, Vector1([1]))
        if np.isnan(dst):
            dst = 10000
        prot_name = "popul" if is_best is None else is_best
        pose.pdb_info().name(prot_name + "_pose_" + str(pdb_id))
        self.pymover.apply(pose)
        rmsd = self.get_rmsd(pose)
        return dst, rmsd, interface, irms

    def render_models(self, pdb_id, SixD_vector, is_best=None, interface=False):
        pose = self.apply_sixD_to_pose(SixD_vector)
        dst = self.scfxn_rosetta.score(pose)
        interface = calc_interaction_energy(pose, self.scfxn_rosetta, Vector1([1]))
        irms = calc_Irmsd(self.native_pose, pose, self.scfxn_rosetta, Vector1([1]))
        if np.isnan(dst):
            dst = 10000
        prot_name = "popul" if is_best is None else is_best
        pose.pdb_info().name(prot_name + "_pose_" + str(pdb_id))
        self.pymover.apply(pose)
        rmsd = self.get_rmsd(pose)
        return dst, rmsd, interface, irms

    def score(self, SixD_vector):
        pose = self.apply_sixD_to_pose(SixD_vector)
        try:
            dst = self.scfxn_rosetta.score(pose)
        except ValueError:
            dst = 10000
        if np.isnan(dst):
            dst = 10000
        return dst

    def size(self):
        return 6

    def convert_positions_to_genotype(self, positions):
        return self.converter.convert_positions_to_genotype(positions)

    def convert_genotype(self, genotype):
        return self.converter.convert_genotype(genotype)

    def get_sol_string(self, sol):
        return " , ".join(["{:.2f}".format(e) for e in sol])

    def apply_sixD_to_pose(self, genotype):
        SixD_vector = self.convert_genotype(genotype)
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
