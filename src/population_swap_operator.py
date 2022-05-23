#!/usr/bin/env python
# coding: utf-8


import random
from pyrosetta import Vector1

from src.utils import get_position_info
from src.individual import Individual
from src.flexbb_swap_operator import FlexbbSwapOperator

from pyrosetta.rosetta.protocols.docking import calc_interaction_energy, calc_Irmsd
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from src.utils import IP_ADDRESS

# SWAP_PROBABILITY = 1.0 # changed from 0.1


class FlexbbSwapOperatorBuilder:
    def __init__(self, config, scfxn):
        self.config = config
        self.scfxn = scfxn

    def build(self):
        if self.config.docking_type_option == "Flexbb":
            self.flexbb_swap_operator = PopulationSwapOperator(self.config, self.scfxn)
        else:
            self.flexbb_swap_operator = None
        return self.flexbb_swap_operator


class PopulationSwapOperator:
    def __init__(self, config, scfxn):
        self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.scfxn = scfxn
        self.scfxn_rosetta = self.scfxn.scfxn_rosetta
        self.config = config
        self.local_search = scfxn.local_search
        self.jobid = self.config.out_path
        self.swap_prob = config.swap_prob

        self.swap_operator = FlexbbSwapOperator(config, scfxn, self.local_search)
        self.log_relax_success = self.jobid + "/relaxed_sucess.csv"
        with open(self.log_relax_success, "w") as file_object:
            file_object.write("#{}\n".format(self.jobid))
        self.log_population_swap = self.jobid + "/population_swap.csv"
        with open(self.log_population_swap, "w") as file_object:
            file_object.write("#{}\n".format(self.jobid))

    def apply(self, population):
        ratio_swap_success = 0
        relaxed_sucess = 0
        new_population = []
        for ind in population:
            if random.uniform(0, 1) <= self.swap_prob:
                pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
                join_pose, idx_r, idx_l, idx_s, relaxed = self.swap_operator.apply_bb_strategy(
                    ind, pose
                )
                self.local_search.local_search_strategy.slide_into_contact.apply(
                    join_pose
                )
                self.local_search.local_search_strategy.docking.apply(join_pose)
                after = self.local_search.energy_score(join_pose)
                if after < ind.score:
                    positions = get_position_info(pose)
                    genotype = self.scfxn.convert_positions_to_genotype(positions)
                    rmsd = CA_rmsd(self.scfxn.native_pose, join_pose)
                    interface = calc_interaction_energy(
                        join_pose, self.scfxn_rosetta, Vector1([1])
                    )
                    irms = calc_Irmsd(
                        self.scfxn.native_pose,
                        join_pose,
                        self.scfxn.scfxn_rosetta,
                        Vector1([1]),
                    )
                    result_individual = Individual(
                        genotype=genotype,
                        score=after,
                        idx_ligand=idx_l,
                        idx_receptor=idx_r,
                        idx_subunit=idx_s,
                        rmsd=rmsd,
                        i_sc=interface,
                        irms=irms,
                    )
                    new_population.append(result_individual)
                    ratio_swap_success += 1
                    if relaxed:
                        relaxed_sucess += 1

                    # print("IMPROVED!!!")
                    if after < self.local_search.best_score:
                        self.local_search.best_score = after
                        self.local_search.best_pose.assign(pose)

                else:
                    new_population.append(ind)
            else:
                new_population.append(ind)
        ratio_swap_success = (
            0 if ratio_swap_success == 0 else ratio_swap_success / len(population)
        )
        with open(self.log_population_swap, "a") as file_object:
            file_object.write("{}\n".format(ratio_swap_success))
        with open(self.log_relax_success, "a") as file_object:
            file_object.write("{}\n".format(relaxed_sucess))

        return self.local_search.best_pose, new_population
