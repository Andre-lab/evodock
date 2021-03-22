import logging
import time

from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.core.scoring import ScoreFunction, ScoreType
from pyrosetta.rosetta.protocols.docking import (DockingLowRes,
                                                 DockingProtocol,
                                                 DockingSlideIntoContact,
                                                 DockMCMCycle, DockMCMProtocol,
                                                 calc_interaction_energy,
                                                 calc_Irmsd)

from src.differential_evolution import Individual
from src.rigid_body_min import CustomPacker, CustomRotamer
from src.utils import get_position_info


def change_axis(initial_position, axes, amount):
    if axes == "x":
        new_xaxis_value = initial_position[3] - amount
        initial_position = (
            initial_position[:3] + [new_xaxis_value] + initial_position[4:]
        )
    if axes == "y":
        new_yaxis_value = initial_position[4] + amount
        initial_position = (
            initial_position[:4] + [new_yaxis_value] + initial_position[5:]
        )
    if axes == "z":
        new_zaxis_value = initial_position[5] + amount
        initial_position = initial_position[:5] + [new_zaxis_value]

    return initial_position


def init_fa_penalization(scfxn):
    init_pose = Pose()
    init_pose.assign(scfxn.dock_pose)
    init_pose.pdb_info().name("init_pose")
    initial_position = get_position_info(init_pose)
    scorefxn_low = ScoreFunction()
    scorefxn_low.set_weight(ScoreType.fa_rep, 0.55)
    amount = 500
    axis = "x"
    positions = change_axis(initial_position, axis, amount)
    genotype = scfxn.convert_positions_to_genotype(positions)
    ind_pose = scfxn.apply_sixD_to_pose(genotype)
    fa_limit = scorefxn_low.score(ind_pose)
    return fa_limit


class LocalSearchPopulation:
    # Options:
    # None: only score and return the poses
    # only_slide: just slide_into_contact
    # custom_rotamer: RotamerTrialsMover
    # custom_packer: stochastic default packer from rosetta
    # mcm_rosetta: mcm protocol mover (high res) from rosetta (2 cycles)
    def __init__(self, scfxn, packer_option="default_combination"):
        self.packer_option = packer_option
        self.scfxn = scfxn
        self.local_logger = logging.getLogger("evodock.local")
        self.local_logger.setLevel(logging.INFO)
        self.slide_into_contact = DockingSlideIntoContact(1)
        if packer_option == "custom_rotamer":
            self.docking = CustomRotamer(scfxn.dock_pose, scfxn.scfxn_rosetta)
        if packer_option == "custom_packer":
            self.docking = CustomPacker(scfxn.dock_pose, scfxn.scfxn_rosetta)
        if packer_option == "custom_pack_rotamer":
            self.docking = CustomPacker(scfxn.dock_pose, scfxn.scfxn_rosetta)
        if packer_option == "mcm_rosetta" or packer_option == "mcm_rosetta_rt":
            mcm_docking = DockMCMProtocol()
            mcm_docking.set_native_pose(scfxn.dock_pose)
            mcm_docking.set_scorefxn(scfxn.scfxn_rosetta)
            if packer_option == "mcm_rosetta_rt":
                mcm_docking.set_rt_min(True)
            else:
                mcm_docking.set_rt_min(False)
            mcm_docking.set_sc_min(False)
            mock_pose = Pose()
            mock_pose.assign(scfxn.dock_pose)
            mcm_docking.apply(mock_pose)
            self.docking = mcm_docking
            self.docking.set_task_factory(mcm_docking.task_factory())
            self.docking.set_ignore_default_task(True)

        if packer_option == "low_res":
            docking = DockingProtocol()
            docking.set_movable_jumps(Vector1([1]))  # set the jump to jump 1
            docking.set_lowres_scorefxn(scfxn.scfxn_rosetta)
            self.docking = docking
            # self.docking = DockingLowRes(scfxn.scfxn_rosetta, Vector1([1]))

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def apply(self, popul):
        start = time.time()
        trials = []
        after_scores = []
        before_scores = []
        if self.packer_option == "None":
            for ind in popul:
                local_search = False
                trial, before, after = self.process_individual(ind, local_search)
                trials.append(trial)
                before_scores.append(before)
                after_scores.append(after)
            self.local_logger.info(
                "avg_without_ls {}".format(sum(before_scores) / len(before_scores))
            )
        else:
            for ind in popul:
                trial, before, after = self.process_individual(ind)
                trials.append(trial)
                before_scores.append(before)
                after_scores.append(after)
            self.local_logger.info(
                "avg_before {}".format(sum(before_scores) / len(before_scores))
            )
            self.local_logger.info(
                "avg_after {}".format(sum(after_scores) / len(after_scores))
            )
        end = time.time()
        self.local_logger.info(" LS time: " + str(end - start))
        return trials

    def set_new_max_translations(self, popul):
        # 1) get all positions
        # 2) set new max_translations
        # 3) convert genotypes to new translation
        all_x = []
        all_y = []
        all_z = []
        all_positions = []
        for ind in popul:
            pose = self.scfxn.apply_sixD_to_pose(ind.genotype)
            positions = get_position_info(pose)
            all_x.append(positions[3])
            all_y.append(positions[4])
            all_z.append(positions[5])
            all_positions.append(positions)

        max_x = max([val for val in all_x])
        max_x += max_x * 0.1
        max_y = max([val for val in all_y])
        max_y += max_y * 0.1
        max_z = max([val for val in all_z])
        max_z += max_z * 0.1

        min_x = min([val for val in all_x])
        min_x += min_x * 0.1
        min_y = min([val for val in all_y])
        min_y += min_y * 0.1
        min_z = min([val for val in all_z])
        min_z += min_z * 0.1

        # print("prev max_trans ")
        # print(self.scfxn.converter.max_trans)
        self.scfxn.converter.max_trans = [max_x, max_y, max_z]
        self.scfxn.converter.min_trans = [min_x, min_y, min_z]
        self.scfxn.converter.bounds = self.scfxn.converter.define_bounds()
        # print("updated max_trans ")
        # print(self.scfxn.converter.max_trans)
        # print(self.scfxn.converter.min_trans)

        for i, ind in enumerate(popul):
            # print("=== before ===")
            # print(all_positions[i])
            ind.genotype = self.scfxn.convert_positions_to_genotype(all_positions[i])
            # print("==== after ===")
            # print(self.scfxn.convert_genotype(ind.genotype))

        return popul

    def process_individual(self, ind, local_search=True):
        pose = self.scfxn.apply_sixD_to_pose(ind)
        before = self.energy_score(pose)
        positions = get_position_info(pose)
        # print(positions)
        if local_search and self.packer_option != "None":
            if self.packer_option != "low_res":
                self.slide_into_contact.apply(pose)
            if self.packer_option != "only_slide":
                self.docking.apply(pose)
            after = self.energy_score(pose)
        else:
            after = before

        rmsd = self.scfxn.get_rmsd(pose)

        interface = calc_interaction_energy(
            pose, self.scfxn.scfxn_rosetta, Vector1([1])
        )
        irms = calc_Irmsd(
            self.scfxn.native_pose, pose, self.scfxn.scfxn_rosetta, Vector1([1])
        )

        # get position from pose
        positions = get_position_info(pose)
        # replace trial with this new positions
        genotype = self.scfxn.convert_positions_to_genotype(positions)
        result_individual = Individual(genotype, after, rmsd, interface, irms)
        # print("THE REAL END")
        return result_individual, before, after
