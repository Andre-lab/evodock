import glob
import logging
import random
import time


from src.individual import Individual
from src.mpi_utils import IndividualMPI
from src.population_swap_operator import FlexbbSwapOperatorBuilder
from src.selection import GreedySelection
from src.single_process import SingleProcessPopulCalculator
from src.mutation_strategies import MutationStrategyBuilder
from src.population import ScorePopulation
from src.trial_generator import (
    TrialGenerator,
    CodeGenerator,
    TriangularGenerator,
    FlexbbTrialGenerator,
)
from src.initialize_population import InitializePopulationBuilder
from src.utils import make_trial
from src.scfxn_fullatom import FAFitnessFunction

# --- MAIN ---------------------------------------------------------------------+


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
        self.max_translation = config.get_max_translation()
        self.scfxn = fitness_function
        self.mutation_strategy = MutationStrategyBuilder(config).build()
        self.flexbb_swap_operator = FlexbbSwapOperatorBuilder(
            config, fitness_function
        ).build()
        # trial_score_config = config
        # trial_score_config.docking_type_option = "Global"
        # self.scfxn = FAFitnessFunction(
        #     self.scfxn.input_pose, self.scfxn.native_pose, trial_score_config
        # )

        self.popul_calculator = ScorePopulation(config, fitness_function)
        self.init_file()

    def init_population(self):
        return InitializePopulationBuilder().run(self)

    def init_file(self):
        header = f"# CONF: maxiter : {self.maxiter}, np : {self.popsize}, f {self.mutate}, cr {self.recombination}\n"
        with open(self.job_id, "w") as file_object:
            file_object.write(header)
            file_object.write("# INIT EVOLUTION\n")
            file_object.write("gen,avg,best,rmsd_from_best\n")
        with open(self.file_time_name, "w") as file_time:
            file_time.write(header)
            file_time.write("# INIT EVOLUTION\n")
            file_time.write("generation_seconds\n")

    def main(self, population):
        self.logger.info(" DE")
        self.population = population
        self.popul_calculator.print_information(self.population)
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

            # --- SCORE KEEPING --------------------------------+
            self.score_keeping()

        return self.best_pdb

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
        self.popul_calculator.print_information(trial_inds, True)
        # self.popul_calculator.pymol_visualization(trial_inds)

        self.previous_gen_scores = self.gen_scores
        (
            self.population,
            self.gen_scores,
            self.trial_scores,
        ) = GreedySelection().apply(trial_inds, self.population)

    def score_keeping(self):
        self.gen_avg = sum(self.gen_scores) / self.popsize
        self.gen_best = min(self.gen_scores)
        self.trial_avg = sum(self.trial_scores) / self.popsize
        self.trial_best = min(self.trial_scores)
        gen_sol = self.population[self.gen_scores.index(min(self.gen_scores))]
        self.logger.info(f"   > GENERATION AVERAGE: {self.gen_avg:.2f}")
        self.logger.info(f"   > GENERATION BEST: {self.gen_best:.2f}")
        self.logger.info(f"   > TRIAL INFO: {self.trial_best:.2f} {self.trial_avg:.2f}")

        (
            self.best_pdb,
            best_SixD_vector,
            best_rmsd,
        ) = self.popul_calculator.render_best(self.generation, gen_sol, self.population)
        best_sol_str = self.popul_calculator.get_sol_string(best_SixD_vector)
        self.logger.info(f"   > BEST SOL: {best_sol_str}")
        self.popul_calculator.print_information(self.population)

        # self.popul_calculator.pymol_visualization(self.population)

        if self.config.out_pdb:
            name = self.config.out_path + "/evolved.pdb"
            self.best_pdb.dump_pdb(name)

        file_object = open(self.job_id, "a")
        evolution_str = f"{self.generation:.0f},{self.gen_avg:.2f},{self.gen_best:.2f},{best_rmsd:.2f}"
        file_object.write(f"{evolution_str}\n")
        file_object.close()
        gen_end = time.time()

        file_time = open(self.file_time_name, "a")
        self.logger.info(f"selection stage in {gen_end - self.gen_start:.2f} s.")
        file_time.write(f"{gen_end - self.gen_start:.2f}\n")
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

    def main(self, population):
        self.logger.info(" DE")
        self.population = population
        self.popul_calculator.print_information(self.population)
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

            # --- Apply Flexbb search ---- #
            self.best_pdb, self.population = self.flexbb_swap_operator.apply(
                self.population
            )

            # --- SCORE KEEPING --------------------------------+
            self.score_keeping()

        return self.best_pdb


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

    def main(self, population):
        self.logger.info(" DE")
        self.population = population
        self.popul_calculator.print_information(self.population)
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

        return self.best_pdb

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
