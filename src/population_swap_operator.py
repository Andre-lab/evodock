#!/usr/bin/env python
# coding: utf-8
import random
from src.utils import get_position_info
from src.individual import Individual
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from cubicsym.utilities import add_id_to_pose_w_base

class FlexbbSwapOperatorBuilder:
    def __init__(self, config, scfxn, dockmetric):
        self.config = config
        self.scfxn = scfxn
        self.dockmetric = dockmetric

    def build(self):
        if self.config.flexbb:
            self.flexbb_swap_operator = PopulationSwapOperator(self.config, self.scfxn, self.dockmetric)
        else:
            self.flexbb_swap_operator = None
        return self.flexbb_swap_operator

class PopulationSwapOperator:
    def __init__(self, config, scfxn, dockmetric):
        self.scfxn = scfxn
        self.scfxn_rosetta = self.scfxn.scfxn_rosetta
        self.dockmetric = dockmetric
        self.config = config
        self.local_search = scfxn.local_search
        self.jobid = self.config.out_path
        self.swap_prob = config.swap_prob
        self.swap_operator = scfxn.local_search.local_search_strategy.swap_operator
        self.log_relax_success = self.jobid + "/relaxed_sucess.csv"
        with open(self.log_relax_success, "w") as file_object:
            file_object.write("gen,#{}\n".format(self.jobid))
        self.log_population_swap = self.jobid + "/population_swap.csv"
        with open(self.log_population_swap, "w") as file_object:
            file_object.write("gen,#{}\n".format(self.jobid))

    def insert_bb_name(self, pose, ind):
        """Retrieve the name of current pose/index in the subunitlist.
        Usefull for keeping track of which current bb are being used by the individual."""
        if self.config.low_memory_mode:
            if is_symmetric(pose):
                ind.subunit_name = self.swap_operator.list_subunits[ind.idx_subunit]
        else:
            raise NotImplementedError("NOT IMPLEMENTED!")

    # def dssp_diff_ok(self, pose, previous_pose):
    #     assert pose.size() == previous_pose.size()
    #     counts = 0
    #     for orig, new in zip(Dssp(previous_pose).get_dssp_secstruct(), Dssp(pose).get_dssp_secstruct()):
    #         if orig != new:
    #             counts += 1
    #     if counts / pose.size() * 100 > self.config.dssp_acceptance:
    #         return False
    #     else:
    #         return True
    #
    # def rmsd_ok(self, pose, previous_pose):
    #     return self.dockmetric.ca_rmsd(pose, previous_pose) <= self.config.rmsd_acceptance
    #
    # def accept_bb_into_ensemble(self, pose, idx_l, idx_r, idx_s, ind):
    #     if self.config.save_bbs:
    #         if is_symmetric(pose):
    #             asym_pose = Pose()
    #             extract_asymmetric_unit(pose, asym_pose, False)
    #             # TODO: check if this still happened
    #             #  Bug in Rosetta where some (atleast 1 I have seen VRT is kept in the asymmetric pose)
    #             if "X" in asym_pose.sequence()m
    #                 print("DEBUG: X FOUND IN ASYMMETRIC POSE!")
    #                 nvrt = asym_pose.sequence().count("X")
    #                 start = asym_pose.size() - nvrt + 1
    #                 assert asym_pose.sequence()[start:] == "X" * nvrt, "The virtual residues are not at the end"
    #                 DeleteRegionMover(start, asym_pose.size()).apply(asym_pose)
    #             previous_pose = self.swap_operator.list_subunits[ind.idx_subunit]
    #             if self.config.low_memory_mode:
    #                 previous_pose = pose_from_file(previous_pose)
    #         else:
    #             raise NotImplementedError
    #         if self.dssp_diff_ok(asym_pose, previous_pose) and self.rmsd_ok(asym_pose, previous_pose):
    #             self.swap_operator.save_bb(pose, idx_l, idx_r, idx_s)

    def apply(self, population, generation):
        ratio_swap_success = 0
        relaxed_sucess = 0
        new_population = []
        for ind in population:
            if random.uniform(0, 1) <= self.swap_prob:
                pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
                join_pose, idx_r, idx_l, idx_s, relaxed = self.swap_operator.apply_bb_strategy(ind, pose)
                add_id_to_pose_w_base(join_pose, idx_s)
                self.local_search.local_search_strategy.slide_into_contact.apply(join_pose)
                self.local_search.local_search_strategy.docking.apply(join_pose)
                after = self.local_search.energy_score(join_pose)
                if after < ind.score:
                    positions = get_position_info(join_pose, self.config.syminfo)
                    genotype = self.scfxn.convert_positions_to_genotype(positions)
                    rmsd = self.dockmetric.ca_rmsd(join_pose)
                    interface = self.dockmetric.interaction_energy(join_pose)
                    irms = self.dockmetric.interface_rmsd(join_pose)
                    result_individual = Individual(
                        idx=ind.idx,
                        genotype=genotype,
                        score=after,
                        idx_ligand=idx_l,
                        idx_receptor=idx_r,
                        idx_subunit=idx_s,
                        rmsd=rmsd,
                        i_sc=interface,
                        irms=irms,
                        flipped=ind.flipped,
                        fixed=ind.fixed
                    )
                    self.insert_bb_name(join_pose, result_individual)
                    new_population.append(result_individual)
                    ratio_swap_success += 1
                    if relaxed:
                        relaxed_sucess += 1
                    if after < self.local_search.best_score:
                        self.local_search.best_score = after
                        self.local_search.best_pose.assign(join_pose)
                else:
                    new_population.append(ind)
            else:
                new_population.append(ind)
        ratio_swap_success = (
            0 if ratio_swap_success == 0 else ratio_swap_success / len(population)
        )
        with open(self.log_population_swap, "a") as file_object:
            file_object.write(f"{generation}," + "{}\n".format(ratio_swap_success))
        with open(self.log_relax_success, "a") as file_object:
            file_object.write(f"{generation}," + "{}\n".format(relaxed_sucess))

        return self.local_search.best_pose, new_population
