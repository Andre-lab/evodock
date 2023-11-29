import logging
import random
import time

import numpy as np

from src.population_swap_operator import FlexbbSwapOperatorBuilder
from src.selection import GreedySelection
from src.mutation_strategies import MutationStrategyBuilder
from src.population import ScorePopulation
from src.trial_generator import (
    TrialGenerator,
    CodeGenerator,
    TriangularGenerator,
    FlexbbTrialGenerator,
)
from src.initialize_population import InitializePopulation
from src.utils import make_trial, get_position_info
from src.landscape_metrics import fitness_distance_correlation
from copy import deepcopy
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from cubicsym.cubicsetup import CubicSetup
from cubicsym.assembly.cubicassembly import CubicSymmetricAssembly

class DifferentialEvolutionAlgorithm:
    def __init__(self, config, fitness_function):
        self.config = config
        self.scheme = config.scheme
        self.logger = logging.getLogger("evodock.de")
        self.logger.setLevel(logging.INFO)
        self.popsize = config.popsize
        self.mutate = config.mutate
        self.recombination = config.recombination
        self.maxiter = config.maxiter
        self.job_id = config.out_path + "/evolution.csv"
        self.ind_size = fitness_function.size()
        self.bounds = [(-1, 1)] * self.ind_size
        self.file_time_name = self.config.out_path + "/time.csv"
        self.max_translation = config.get_max_translation(fitness_function.initial_pose)
        self.scfxn = fitness_function
        self.mutation_strategy = MutationStrategyBuilder(config).build(self.scfxn.size())
        self.flexbb_swap_operator = FlexbbSwapOperatorBuilder(
            config, fitness_function, self.scfxn.dockmetric,
        ).build()
        if fitness_function.native_pose is not None:
            if self.config.syminfo:
                self.optimal_solution = get_position_info(fitness_function.native_symmetric_pose, self.config.syminfo)
            else:
                self.optimal_solution = get_position_info(fitness_function.native_pose, self.config.syminfo)
        else:
            self.optimal_solution = None
        self.popul_calculator = ScorePopulation(config, fitness_function)
        self.init_file()
        self.all_docks = {"ind": [], "genotype": [], "i_sc": [], "score": []}

    def init_population(self):
        self.population = InitializePopulation(self.config, self.logger, self.popul_calculator, self.scfxn).init_population()

    def init_file(self):
        # header = f"# CONF: maxiter : {self.maxiter}, np : {self.popsize}, f {self.mutate}, cr {self.recombination}\n"
        with open(self.job_id, "w") as file_object:
            file_object.write("gen,avg,best,rmsd_from_best,fdc\n")
        with open(self.file_time_name, "w") as file_time:
            file_time.write("gen,generation_seconds\n")

    def main(self):
        self.logger.info(" DE")
        self.popul_calculator.print_information(self.population, "init")
        self.popul_calculator.print_genotype_values(self.population, "init")
        for generation in range(1, self.maxiter + 1):
            self.generation = generation
            self.logger.info(f" GENERATION: {self.generation}")
            self.gen_start = time.time()
            if self.config.selection == "total":
                self.gen_scores = [ind.score for ind in self.population]
            else:
                self.gen_scores = [ind.i_sc for ind in self.population]
            # cycle through each individual in the population

            # --- TRIAL CREATION (step #3.A) -------------+
            trials = self.generate_trial_population()

            # --- SELECTION (step #3.C) -------------+
            self.selection(trials)

            # --- SCORE KEEPING --------------------------------+
            self.score_keeping()

    def generate_trial_population(self):
        trials = []
        trial_generator = TrialGenerator(
            self.config, self.popul_calculator.local_search
        )
        for j in range(0, self.popsize):
            # --- MUTATION (step #3.A) ---------------------+
            # select 3 random vector index positions [0, self.popsize)
            # not including current vector (j)
            v_trial = trial_generator.build(j, self.population, self.gen_scores)
            trials.append(v_trial)
        return trials

    def selection(self, trials):
        trial_inds = trials
        self.popul_calculator.print_information(trial_inds, self.generation, True)
        # self.popul_calculator.pymol_visualization(trial_inds)

        self.previous_gen_scores = self.gen_scores
        (
            self.population,
            self.gen_scores,
            self.trial_scores,
        ) = GreedySelection().apply(trial_inds, self.population, self.config.selection)

    def map_all_docks(self, population):
        for ind in population:
            self.all_docks["ind"].append(deepcopy(ind))
            self.all_docks["genotype"].append(deepcopy(ind.genotype))
            self.all_docks["i_sc"].append(deepcopy(ind.i_sc))
            self.all_docks["score"].append(deepcopy(ind.score))

    def score_keeping(self):
        self.gen_avg = sum(self.gen_scores) / self.popsize
        self.gen_best = min(self.gen_scores)
        self.trial_avg = sum(self.trial_scores) / self.popsize
        self.trial_best = min(self.trial_scores)
        best_ind_solution = self.population[self.gen_scores.index(min(self.gen_scores))]
        self.logger.info(f"   > GENERATION AVERAGE: {self.gen_avg:.2f}")
        self.logger.info(f"   > GENERATION BEST: {self.gen_best:.2f}")
        self.logger.info(f"   > TRIAL INFO: {self.trial_best:.2f} {self.trial_avg:.2f}")

        best_SixD_vector, best_rmsd = self.popul_calculator.render_best(self.generation, best_ind_solution, self.population)

        self.map_all_docks(self.population)

        self.popul_calculator.print_genotype_values(self.population, self.generation)

        best_sol_str = self.popul_calculator.get_sol_string(best_SixD_vector)
        improved = np.array(
            [y - x for x, y in zip(self.previous_gen_scores, self.gen_scores) if y < x]
        )

        avg_improved = np.average(improved) if len(improved) > 0 else 0
        self.logger.info(
            f" improved {len(improved)} individuals with average {avg_improved:.2f}"
        )

        self.logger.info(f"   > BEST SOLUTION: {best_sol_str}")
        self.popul_calculator.print_information(self.population, self.generation)

        # print ensemble names
        self.popul_calculator.print_ensembles(self.population, self.generation)

        # self.popul_calculator.pymol_visualization(self.population)

        if self.config.out_pdb or self.config.output_pdb_per_generation:
            if self.config.out_pdb:
                name = self.config.out_path + "/evolved.pdb"
            elif self.config.output_pdb_per_generation:
                name = self.config.out_path + f"/evolved_{self.generation}.pdb"
            if self.config.flexbb:
                best_pose, _, _, _ = self.scfxn.local_search.local_search_strategy.get_pose(best_ind_solution)
            else:
                best_pose = self.scfxn.local_search.local_search_strategy.get_pose(best_ind_solution)
            if is_symmetric(best_pose):
                cs = CubicSetup(self.config.syminfo.input_symdef)
                if cs.is_cubic():
                    # output full symmetric structure
                    CubicSymmetricAssembly.from_pose_input(best_pose, cs).output(filename=name, format="cif")
            else:
                best_pose.dump_pdb(name)

        if self.optimal_solution is not None:
            fdc = fitness_distance_correlation(
                self.population, self.optimal_solution, self.scfxn
            )
        else:
            fdc = -1

        file_object = open(self.job_id, "a")
        evolution_str = f"{self.generation:.0f},{self.gen_avg:.2f},{self.gen_best:.2f},{best_rmsd:.2f},{fdc:.2f}"
        file_object.write(f"{evolution_str}\n")
        file_object.close()
        gen_end = time.time()

        file_time = open(self.file_time_name, "a")
        self.logger.info(f"selection stage in {gen_end - self.gen_start:.2f} s.")
        file_time.write(f"{self.generation},{gen_end - self.gen_start:.2f}\n")
        file_time.close()

    def apply_popul_flexbb(self, population):
        return self.flexbb_swap_operator.apply(population)


class FlexbbDifferentialEvolution(DifferentialEvolutionAlgorithm):


    def generate_trial_population(self):
        trials = []
        trial_generator = FlexbbTrialGenerator(
            self.config, self.popul_calculator.local_search
        )
        for j in range(0, self.popsize):
            # --- MUTATION (step #3.A) ---------------------+
            # select 3 random vector index positions [0, self.popsize)
            # not including current vector (j)
            v_trial = trial_generator.build(j, self.population, self.gen_scores)
            trials.append(v_trial)
        return trials

    def main(self):
        self.logger.info(" DE")
        self.popul_calculator.print_information(self.population, gen="init")
        self.popul_calculator.print_ensembles(self.population, gen="init")
        self.popul_calculator.print_flip_fix_info(self.population)
        self.popul_calculator.print_genotype_values(self.population, gen="init")
        self.archive_restart = [0] * self.popsize
        for generation in range(1, self.maxiter + 1):
            self.generation = generation
            self.logger.info(f" GENERATION: {self.generation}")
            self.gen_start = time.time()
            if self.config.selection == "total":
                self.gen_scores = [ind.score for ind in self.population]
            else:
                self.gen_scores = [ind.i_sc for ind in self.population]
            # cycle through each individual in the population

            # --- TRIAL CREATION (step #3.A) -------------+
            trials = self.generate_trial_population()

            # --- SELECTION (step #3.C) -------------+
            self.selection(trials)

            # --- Apply Flexbb search ---- #
            self.population = self.flexbb_swap_operator.apply(self.population, self.generation)

            # --- SCORE KEEPING --------------------------------+
            self.score_keeping()


class TriangularDE(DifferentialEvolutionAlgorithm):
    def generate_trial_population(self):
        rmax = 1
        rmin = 0.1
        rgen = rmax - ((self.generation / self.maxiter) * (rmax - rmin))
        trials = []
        trial_generator = TriangularGenerator(
            self.config, self.generation, self.popul_calculator.local_search
        )
        for j in range(0, self.popsize):
            # --- MUTATION (step #3.A) ---------------------+
            # select 3 random vector index positions [0, self.popsize)
            # not including current vector (j)
            v_trial = trial_generator.build(j, self.population, self.gen_scores)
            trials.append(v_trial)
        return trials

    def main(self):
        self.logger.info(" DE")
        self.popul_calculator.print_information(self.population, None)
        self.archive_restart = [0] * self.popsize
        for generation in range(1, self.maxiter + 1):
            self.generation = generation
            self.logger.info(f" GENERATION: {self.generation}")
            self.gen_start = time.time()
            self.gen_scores = [ind.score for ind in self.population]
            # cycle through each individual in the population

            # --- TRIAL CREATION (step #3.A) -------------+
            trials = self.generate_trial_population()

            # --- SELECTION (step #3.C) -------------+
            self.selection(trials)

            # --- RESTART ---- #
            self.restart_mechanism()

            # --- SCORE KEEPING --------------------------------+
            self.score_keeping()

    def restart_mechanism(self):
        self.gen_best = min(self.gen_scores)
        for j in range(0, self.popsize):
            if abs(self.gen_scores[j] - self.previous_gen_scores[j]) < 0.0001:
                self.archive_restart[j] += 1
            else:
                self.archive_restart[j] = 0
            if self.archive_restart[j] == 15:
                if self.gen_scores[j] != self.gen_best:
                    self.archive_restart[j] = 0
                    indv = self.population[j].genotype
                    jrand = random.randrange(0, len(self.bounds))
                    indv[jrand] = random.uniform(
                        self.bounds[jrand][0], self.bounds[jrand][1]
                    )
                    ind = make_trial(j, indv, 0, 0)
                    (
                        ind,
                        before,
                        after,
                    ) = self.scfxn.local_search.process_individual(ind)
                    self.population[j] = ind
                    self.gen_scores[j] = after
