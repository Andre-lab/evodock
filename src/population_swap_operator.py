#!/usr/bin/env python
# coding: utf-8



import random
from pyrosetta import Pose, Vector1

from src.utils import get_position_info
from src.individual import Individual
from src.flexbb_swap_operator import FlexbbSwapOperator
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.pose import remove_virtual_residues

from pyrosetta.rosetta.protocols.docking import calc_interaction_energy, calc_Irmsd
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from src.utils import IP_ADDRESS

SWAP_PROBABILITY = 0.1

            
class PopulationSwapOperator:
    def __init__(self, config, scfxn, local_search):
        self.pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        self.scfxn = scfxn
        self.scfxn_rosetta = self.scfxn.scfxn_rosetta
        self.config = config
        self.local_search = local_search
        self.jobid = self.config.jobid

        self.swap_operator = FlexbbSwapOperator(config, scfxn, local_search)
        self.log_relax_success = self.jobid.replace("evolution", "relaxed_sucess")
        with open(self.log_relax_success, "w") as file_object:
            file_object.write("#{}\n".format(self.jobid))
        self.log_population_swap = self.jobid.replace("evolution", "population_swap")
        with open(self.log_population_swap, "w") as file_object:
            file_object.write("#{}\n".format(self.jobid))
                

    def apply(self, population):
        ratio_swap_success = 0
        relaxed_sucess = 0
        new_population = []
        for ind in population:
            if random.uniform(0, 1) < SWAP_PROBABILITY:
                pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
                join_pose, idx_r, idx_l, relaxed = self.swap_operator.apply_bb_strategy(ind, pose)
                self.local_search.local_search_strategy.slide_into_contact.apply(join_pose)
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
                        self.scfxn.native_pose, join_pose, self.scfxn.scfxn_rosetta, Vector1([1]),
                    )
                    result_individual = Individual(
                        genotype,
                        after,
                        idx_l,
                        idx_r,
                        rmsd,
                        interface,
                        irms,
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
        ratio_swap_success = 0 if ratio_swap_success == 0 else ratio_swap_success/len(population)
        with open(self.log_population_swap, "a") as file_object:
            file_object.write("{}\n".format(ratio_swap_success))
        with open(self.log_relax_success, "a") as file_object:
            file_object.write("{}\n".format(relaxed_sucess))
        

        
        return self.local_search.best_pose, new_population