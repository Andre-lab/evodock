#!usr/bin/env python

import logging
from math import sqrt

import numpy as np
from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreFunctionFactory
from pyrosetta.rosetta.protocols.docking import (calc_interaction_energy,
                                                 calc_Irmsd)
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from scipy.spatial.transform import Rotation as R

from src.genotype_converter import (GlobalGenotypeConverter,
                                    LocalGenotypeConverter)
from src.position_utils import build_axis, to_rosetta
from src.utils import IP_ADDRESS

ref2015_terms = {
    "fa_atr": 1,
    "fa_rep": 0.55,
    "fa_sol": 1.0,
    "fa_intra_sol_xover4": 1.0,
    "lk_ball_wtd": 1.0,
    "fa_intra_rep": 0.005,
    "fa_elec": 1.0,
    "pro_close": 1.25,
    "hbond_sr_bb": 1.0,
    "hbond_lr_bb": 1.0,
    "hbond_bb_sc": 1.0,
    "hbond_sc": 1.0,
    "dslf_fa13": 1.25,
    "rama_prepro": 0.45,
    "omega": 0.4,
    "p_aa_pp": 0.6,
    "fa_dun": 0.7,
    "yhh_planarity": 0.625,
    "ref": 1,
}


def l2_norm(zd_reference, zd_current):
    L2_norm = 0
    for i, v in enumerate(zd_current):
        L2_norm += (zd_current[i] - zd_reference[i]) * (zd_current[i] - zd_reference[i])
    sqrt_L2_norm = sqrt(L2_norm)
    return sqrt_L2_norm


class FAFitnessFunction:
    def __init__(self, native_pose, input_pose, trans_max_magnitude, refinement=False):
        self.logger = logging.getLogger("evodock.scfxn")
        self.native_pose = native_pose
        self.input_pose = input_pose
        self.logger.setLevel(logging.INFO)
        self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.scfxn_rosetta = ScoreFunctionFactory.create_score_function("ref2015")
        self.dock_pose = Pose()
        if refinement:
            self.converter = LocalGenotypeConverter(input_pose)
        else:
            self.converter = GlobalGenotypeConverter(input_pose)

        self.pymover.apply(input_pose)
        self.trans_max_magnitude = trans_max_magnitude
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
            # dst = calc_interaction_energy(pose, self.scfxn_rosetta, Vector1([1]))
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

    def convert_positions_to_genotype(self, positions):
        return self.converter.convert_positions_to_genotype(positions)

    def convert_genotype(self, genotype):
        return self.converter.convert_genotype(genotype)

    def get_dict_scores(self, pose):
        pose_scores = []
        energies = pose.energies().total_energies()
        energies_as_list = [i.strip("( )") for i in str(energies).split(") (")]

        dict_score = {}
        for e in energies_as_list:
            term, unweighted_val = e.split()
            term = term.replace(";", "")
            if term in ref2015_terms:
                weighted_val = float(unweighted_val) * ref2015_terms[term]
                dict_score[term] = weighted_val
                pose_scores.append(": ".join([term, str(weighted_val)]))

        return dict_score

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
