#!usr/bin/env python
import copy
import logging
import numpy as np
import random
from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
# from pyrosetta.rosetta.protocols.moves import PyMOLMover
from scipy.spatial.transform import Rotation as R
from src.genotype_converter import GlobalGenotypeConverter
from src.position_utils import build_axis, to_rosetta
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from src.utils import get_rotation_euler, get_translation
from src.local_search import LocalSearchPopulation
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from src.individual import Individual
from pyrosetta.rosetta.numeric import xyzVector_double_t
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory, ScoreTypeManager
from symmetryhandler.reference_kinematics import set_jumpdof_int_int

class FAFitnessFunction:
    def __init__(self, input_pose, native_pose, config, dockmetric, native_symmetric_pose=None):
        self.logger = logging.getLogger("evodock.scfxn")
        if  is_symmetric(input_pose):
            self.trans_max_magnitude = None
        else:
            self.trans_max_magnitude = config.get_max_translation()
        self.native_pose = native_pose
        self.native_symmetric_pose = native_symmetric_pose
        self.original_reference_scores = []
        self.config = config
        self.dockmetric = dockmetric
        # self.native_pose = Pose()
        # self.native_pose.assign(native_pose)
        # self.input_pose = Pose()
        # self.input_pose.assign(input_pose)
        self.logger.setLevel(logging.INFO)
        # self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.scfxn_rosetta = ScoreFunctionFactory.create_score_function("ref2015")
        self.dock_pose = Pose()
        self.dock_pose.assign(input_pose)
        self.dock_pose.pdb_info().name("INIT_STATE")
        if config.syminfo:
            SetupForSymmetryMover(config.syminfo.input_symdef).apply(self.dock_pose)
            self.config.syminfo.cubicboundary.turn_on_constraint_for_score(self.scfxn_rosetta)

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
        self.slide_attempts = None
        self.setup_interface_score()
        self.local_search = LocalSearchPopulation(
            self, config.local_search_option, config, dockmetric
        )

    def setup_interface_score(self):
        self.in_vdw_score = ScoreFunctionFactory.create_score_function("empty")
        scoretype = ScoreTypeManager().score_type_from_name("interchain_vdw")
        self.in_vdw_score.set_weight(scoretype, 1.0)

    def slide_away(self, pose, rnd_trans_mag=500, max_attemps=50, min_vdw_score=0.05):
        """Apply translations to the pose in order to move all chains out of the way of each other. To make sure the
        chains are not in contact the interchain_vdw score is evaluated. If it is above the acceptable threshold an increase in the
        translational magnitude it applied randomly to any of the translation dofs."""
        direction = 1
        self.slide_attempts = 1
        if self.syminfo:
            while self.slide_attempts <= max_attemps:
                for jumpid, dofid, transmag in self.syminfo.normalize_trans_map:
                    jump = pose.jump(jumpid)
                    trans = get_translation(jump)
                    trans[dofid - 1] += transmag * direction
                    jump.set_translation(xyzVector_double_t(*trans))
                    pose.set_jump(jumpid, jump)
                vdw_score = self.in_vdw_score.score(pose)
                if vdw_score > min_vdw_score:
                    # reset jumps by sliding back
                    self.slide_in(pose)
                    # increase translation in the trans_map
                    trans_map = random.choice(self.syminfo.normalize_trans_map)
                    trans_map[-1] += rnd_trans_mag
                    self.slide_attempts += 1
                else:
                    break
        else:
            raise NotImplementedError

    def slide_in(self, pose):
        """Apply translations to the pose in order to move all chains back to their original position. This functions is the same as
        the slide_away function but without clash detection. """
        direction = -1
        if self.syminfo:
            for jumpid, dofid, transmag in self.syminfo.normalize_trans_map:
                jump = pose.jump(jumpid)
                trans = get_translation(jump)
                trans[dofid - 1] += transmag * direction
                jump.set_translation(xyzVector_double_t(*trans))
                pose.set_jump(jumpid, jump)
        else:
            raise NotImplementedError

    def reset_trans_map(self):
        """resets the normalize_trans_map to its original value"""
        self.syminfo.normalize_trans_map = copy.deepcopy(self.syminfo.org_normalize_trans_map)

    def score(self, pose, ind: Individual):
        """Scores the pose and if normalize_scores is set it will normalize scores accordingly."""
        if self.normalize_score:
            raise NotImplementedError # FIXME: implement properly or delete
            # if self.syminfo:
            #     self.slide_away(pose)
            #     monomeric_score = self.scfxn_rosetta.score(pose)
            #     if monomeric_score > 0:
            #         from pyrosetta.rosetta.protocols.minimization_packing import RotamerTrialsMover, MinMover, PackRotamersMover
            #         from src.symmetry import SymDockMCMCycle
            #         tf = SymDockMCMCycle.create_taskfactory(self.native_pose)
            #         packer = PackRotamersMover(self.scfxn_rosetta)
            #         packer.task_factory(tf)
            #         packer.apply(pose)
            #     self.slide_in(pose)
            #     multimeric_score = self.scfxn_rosetta.score(pose)
            #     self.reset_trans_map()
            #     return multimeric_score - monomeric_score
            # else:
            #     raise NotImplementedError
        else:
            return self.scfxn_rosetta.score(pose)

    def get_rmsd(self, pose):
        rmsd = self.dockmetric.ca_rmsd()
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
            dofsvector = iter(DoFs_vector)
            for jump, dofs in zip(self.syminfo.jumps_int, self.syminfo.dofs_int):
                for dof in dofs:
                    set_jumpdof_int_int(ind_pose, jump, dof, next(dofsvector))
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
        interface = self.dockmetric.interaction_energy(pose)
        irms = self.dockmetric.i_rmsd(pose)
        if np.isnan(dst):
            dst = 10000
        prot_name = "popul" if is_best is None else is_best
        pose.pdb_info().name(prot_name + "_pose_" + str(pdb_id))
        # self.pymover.apply(pose)
        rmsd = self.dockmetric.ca_rmsd(pose)
        return dst, rmsd, interface, irms

    def render_models(self, pdb_id, genotype, is_best=None, interface=False):
        pose = self.apply_genotype_to_pose(genotype)
        dst = self.scfxn_rosetta.score(pose)
        interface = self.dockmetric.interaction_energy(pose)
        irms = self.dockmetric.i_rmsd(pose)
        if np.isnan(dst):
            dst = 10000
        prot_name = "popul" if is_best is None else is_best
        pose.pdb_info().name(prot_name + "_pose_" + str(pdb_id))
        # self.pymover.apply(pose)
        rmsd = self.dockmetric.ca_rmsd(pose)
        return dst, rmsd, interface, irms
