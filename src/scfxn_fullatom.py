#!usr/bin/env python
import logging
import numpy as np
from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreFunctionFactory, CA_rmsd_symmetric
from pyrosetta.rosetta.protocols.docking import (calc_interaction_energy,
                                                 calc_Irmsd)
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreFunctionFactory
from pyrosetta.rosetta.protocols.docking import (
    calc_interaction_energy,
    calc_Irmsd,
)
# from pyrosetta.rosetta.protocols.moves import PyMOLMover
from scipy.spatial.transform import Rotation as R
from src.genotype_converter import GlobalGenotypeConverter
from src.position_utils import build_axis, to_rosetta
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from src.utils import get_rotation_euler, get_translation
from src.local_search import LocalSearchPopulation
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.std import istringstream
from pyrosetta.rosetta.core.conformation.symmetry import SymmData
from src.individual import Individual

class FAFitnessFunction:
    def __init__(self, input_pose, native_pose, config):
        self.logger = logging.getLogger("evodock.scfxn")
        self.trans_max_magnitude = config.get_max_translation()
        self.native_pose = native_pose
        self.original_reference_scores = []
        # self.native_pose = Pose()
        # self.native_pose.assign(native_pose)
        # self.input_pose = Pose()
        # self.input_pose.assign(input_pose)
        self.logger.setLevel(logging.INFO)
        # self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.scfxn_rosetta = ScoreFunctionFactory.create_score_function("ref2015")
        self.dock_pose = Pose()
        self.dock_pose.assign(input_pose)
        if config.syminfo:
            SetupForSymmetryMover(config.syminfo.get("input_symdef")).apply(self.dock_pose)
        self.dock_pose.pdb_info().name("INIT_STATE")

        self.converter = GlobalGenotypeConverter(
            self.dock_pose, self.trans_max_magnitude, config.syminfo
        )
        # self.pymover.apply(self.dock_pose)
        self.jump_num = 1
        self.ax1, self.ax2, self.ax3 = build_axis()
        mros_temp = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.mros, _ = to_rosetta(mros_temp, [0, 0, 0])
        self.syminfo = config.syminfo
        self.normalize_score = config.normalize_score
        self.local_search = LocalSearchPopulation(
            self, config.local_search_option, config
        )

    def normalize_scores(self):
        if self.normalize_score:
            assert self.original_reference_scores, "self.reference_scores is not filled with values yet!"
            if self.syminfo:
                s = SymmData()
                s.read_symmetry_data_from_file(self.syminfo["input_symdef"])
                multiplier = s.get_score_multiply_subunit()[s.get_score_subunit()]
                self.original_reference_scores = [multiplier * i for i in self.original_reference_scores]
                min_val = min(self.original_reference_scores)
                self.reference_scores_subunit = [min_val - i for i in self.original_reference_scores]
            else:
                raise NotImplementedError

    def score(self, pose, ind: Individual):
        """Scores the pose and if normalize_scores is set it will normalize scores accordingly."""
        score = self.scfxn_rosetta.score(pose)
        if self.normalize_score:
            if self.syminfo:
                # IT IS VERY IMPORTANT WE '+' BECAUSE THE REFERENCE VALUES ARE SUPPOSED TO SUBTRACT
                return score + self.reference_scores_subunit[ind.idx_subunit]
            else:
                raise NotImplementedError
        else:
            return score

    def get_rmsd(self, pose):
        rmsd = CA_rmsd(self.native_pose, pose)
        return rmsd

    def size(self):
        return self.converter.size

    def convert_positions_to_genotype(self, positions):
        return self.converter.convert_positions_to_genotype(positions)

    def convert_genotype_to_positions(self, genotype):
        return self.converter.convert_genotype(genotype)

    def get_sol_string(self, sol):
        return " , ".join(["{:.2f}".format(e) for e in sol])

    def get_solution_from_positions(self, DoFs_vector):
        ind_pose = Pose()
        ind_pose.assign(self.dock_pose)
        euler = np.asarray(DoFs_vector[0:3])
        r = R.from_euler("xyz", euler, degrees=True).as_matrix()
        flexible_jump = ind_pose.jump(ind_pose.num_jump())
        rosetta_rotation, rosetta_translation = to_rosetta(r, DoFs_vector[3:])
        flexible_jump.set_rotation(rosetta_rotation)
        flexible_jump.set_translation(rosetta_translation)
        ind_pose.set_jump(ind_pose.num_jump(), flexible_jump)
        # now is time to score the joined pose (ind_pose)
        return ind_pose

    def apply_genotype_to_pose(self, genotype):
        DoFs_vector = self.convert_genotype_to_positions(genotype)
        ind_pose = Pose()
        ind_pose.assign(self.dock_pose)
        if is_symmetric(ind_pose):
            Dofs_vector = iter(DoFs_vector)
            for jump, dofs in zip(self.syminfo.get("jumps_int"), self.syminfo.get("dofs_int")):
                flexible_jump = ind_pose.jump(jump)
                rot = get_rotation_euler(flexible_jump)
                trans = get_translation(flexible_jump)
                for dof in dofs:
                    if dof < 4:
                        trans[dof - 1] = next(Dofs_vector)
                    else:
                        rot[dof - 4] = next(Dofs_vector)
                rot = R.from_euler("xyz", rot, degrees=True).as_matrix()
                rosetta_rotation, rosetta_translation = to_rosetta(rot, trans)
                flexible_jump.set_rotation(rosetta_rotation)
                flexible_jump.set_translation(rosetta_translation)
                ind_pose.set_jump(jump, flexible_jump)
        else:
            # TODO: dont want to add self.jump_num instead of ind_pose.num_jump() ?
            euler = np.asarray(DoFs_vector[0:3])
            r = R.from_euler("xyz", euler, degrees=True).as_matrix()
            flexible_jump = ind_pose.jump(ind_pose.num_jump())
            rosetta_rotation, rosetta_translation = to_rosetta(r, DoFs_vector[3:])
            flexible_jump.set_rotation(rosetta_rotation)
            flexible_jump.set_translation(rosetta_translation)
            ind_pose.set_jump(ind_pose.num_jump(), flexible_jump)
        return ind_pose

    def render_individual(self, pdb_id, individual, is_best=None, interface=False):
        pose = individual.pose
        dst = self.scfxn_rosetta.score(pose)
        interface = calc_interaction_energy(pose, self.scfxn_rosetta, Vector1([1]))
        irms = calc_Irmsd(self.native_pose, pose, self.scfxn_rosetta, Vector1([1]))
        if np.isnan(dst):
            dst = 10000
        prot_name = "popul" if is_best is None else is_best
        pose.pdb_info().name(prot_name + "_pose_" + str(pdb_id))
        # self.pymover.apply(pose)
        rmsd = self.get_rmsd(pose)
        return dst, rmsd, interface, irms

    def render_models(self, pdb_id, genotype, is_best=None, interface=False):
        pose = self.apply_genotype_to_pose(genotype)
        dst = self.scfxn_rosetta.score(pose)
        interface = calc_interaction_energy(pose, self.scfxn_rosetta, Vector1([1]))
        irms = calc_Irmsd(self.native_pose, pose, self.scfxn_rosetta, Vector1([1]))
        if np.isnan(dst):
            dst = 10000
        prot_name = "popul" if is_best is None else is_best
        pose.pdb_info().name(prot_name + "_pose_" + str(pdb_id))
        # self.pymover.apply(pose)
        rmsd = self.get_rmsd(pose)
        return dst, rmsd, interface, irms
