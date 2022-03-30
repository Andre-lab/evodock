import glob
import logging
import random
import time

from src.genotype_converter import RefineCluspro, generate_genotype
from src.individual import Individual
from src.mpi_utils import IndividualMPI
from src.population_swap_operator import PopulationSwapOperator
from src.selection import GreedySelection
from src.single_process import SingleProcessPopulCalculator
from src.mutation_strategies import MutationStrategyBuilder

from src.trial_generator import TrialGenerator, CodeGenerator, TriangularGenerator

# --- MAIN ---------------------------------------------------------------------+


def make_trial(idx, genotype, ligand, receptor):
    ind = IndividualMPI(idx, genotype)
    ind.genotype = genotype
    ind.score = 1000
    ind.idx_ligand = ligand
    ind.idx_receptor = receptor
    ind.rmsd = 0
    ind.i_sc = 0
    ind.irms = 0
    return ind


class DifferentialEvolutionAlgorithm:
    def __init__(self, popul_calculator, config):
        self.config = config
        self.scheme = config.scheme
        self.logger = logging.getLogger("evodock.de")
        self.logger.setLevel(logging.INFO)
        self.popul_calculator = popul_calculator
        self.popsize = config.popsize
        self.mutate = config.mutate
        self.recombination = config.recombination
        self.maxiter = config.maxiter
        self.job_id = config.out_path + "/evolution.log"
        self.ind_size = popul_calculator.cost_func.size()
        self.bounds = [(-1, 1)] * self.ind_size
        self.file_time_name = self.config.out_path + "/time.log"
        self.max_translation = config.get_max_translation()
        self.mutation_strategy = MutationStrategyBuilder(config).build()
        if (
            self.config.docking_type_option == "Unbound"
            and self.config.bb_strategy == "popul_library"
        ):
            self.flexbb_swap_operator = PopulationSwapOperator(
                config, self.popul_calculator.scfxn, self.popul_calculator.local_search
            )
        else:
            self.flexbb_swap_operator = None
        self.init_file()

    def init_population(self, popsize=None, docking_type="Global"):
        # --- INITIALIZE A POPULATION (step #1) ----------------+

        # population_calculator = SingleProcessPopulCalculator(
        #     self.popul_calculator.cost_func, self.config
        # )

        population_calculator = self.popul_calculator

        if docking_type == "RefineCluspro":
            refCluspro = RefineCluspro(self.config, self.max_translation)

        if popsize is None:
            popsize = self.popsize
        self.logger.info(" init population")
        population = []
        popul = []
        for i in range(0, popsize):
            indv = []
            if docking_type == "Global":
                for j in range(len(self.bounds)):
                    indv.append(random.uniform(self.bounds[j][0], self.bounds[j][1]))

            else:
                if docking_type == "RefineCluspro":
                    indv = refCluspro.refine_cluspro(i)
                else:
                    indv = generate_genotype(
                        self.popul_calculator.scfxn.input_pose, self.max_translation
                    )

            # todo: read file and count *.pdb

            if self.config.docking_type_option == "Unbound":
                len_receptors = len(glob.glob(self.config.path_receptors)) - 1
                len_ligands = len(glob.glob(self.config.path_ligands)) - 1
                idx_receptor = random.randint(0, len_receptors)
                idx_ligand = random.randint(0, len_ligands)
            else:
                idx_ligand, idx_receptor = 1, 1

            popul.append(make_trial(i, indv, idx_ligand, idx_receptor))
            population.append(Individual(indv, idx_ligand, idx_receptor, 0, 1000))

        init_population = True
        start = time.time()
        population = population_calculator.run(popul, init_population)
        end = time.time()
        self.logger.info(f"population init in {end - start} seconds")
        # self.popul_calculator.cost_func.pymol_visualization(population)
        return population

    def init_file(self):
        with open(self.job_id, "w") as file_object:
            file_object.write(
                "CONF: maxiter : {}, np : {}, f {}, cr {} \n".format(
                    self.maxiter, self.popsize, self.mutate, self.recombination
                )
            )
            file_object.write("INIT EVOLUTION\n")
        with open(self.file_time_name, "w") as file_time:
            file_time.write(
                "CONF: maxiter : {}, np : {}, f {}, cr {} \n".format(
                    self.maxiter, self.popsize, self.mutate, self.recombination
                )
            )
            file_time.write("INIT EVOLUTION\n")

    def main(self, population):
        self.logger.info(" DE")
        # --- SOLVE --------------------------------------------+

        self.popul_calculator.cost_func.print_information(population)
        # _ = input("continue > ")
        # cycle through each generation (step #2)
        archive_restart = [0] * self.popsize
        for i in range(1, self.maxiter + 1):
            self.logger.info(" GENERATION:" + str(i))
            gen_start = time.time()
            file_object = open(self.job_id, "a")
            file_time = open(self.file_time_name, "a")

            file_object.write("GENERATION: \t" + str(i) + "\t")
            gen_scores = [ind.score for ind in population]

            gen_sol = population[gen_scores.index(min(gen_scores))]
            # cycle through each individual in the population

            rmax = 1
            rmin = 0.1
            rgen = rmax - ((i / self.maxiter) * (rmax - rmin))

            trials = []
            trial_generator = TriangularGenerator(
                self.config, i, self.popul_calculator.cost_func.local_search
            )
            for j in range(0, self.popsize):
                # --- MUTATION (step #3.A) ---------------------+
                # select 3 random vector index positions [0, self.popsize)
                # not including current vector (j)
                v_trial = trial_generator.build(j, population, gen_scores)
                trials.append(v_trial)

            # --- SELECTION (step #3.C) -------------+

            trial_inds = trials
            self.popul_calculator.cost_func.print_information(trial_inds, True)
            # self.popul_calculator.cost_func.pymol_visualization(trial_inds)

            previous_gen_scores = gen_scores
            population, gen_scores, trial_scores = GreedySelection().apply(
                trial_inds, population
            )

            # --- RESTART ---- #
            gen_best = min(gen_scores)
            # for j in range(0, self.popsize):
            if False:
                if abs(gen_scores[j] - previous_gen_scores[j]) < 0.0001:
                    archive_restart[j] += 1
                else:
                    archive_restart[j] = 0
                if archive_restart[j] == 15:
                    if gen_scores[j] != gen_best:
                        archive_restart[j] = 0
                        indv = population[j].genotype
                        jrand = random.randrange(0, len(self.bounds))
                        indv[jrand] = random.uniform(
                            self.bounds[jrand][0], self.bounds[jrand][1]
                        )
                        ind = make_trial(j, indv, 0, 0)
                        (
                            ind,
                            before,
                            after,
                        ) = self.popul_calculator.cost_func.local_search.process_individual(
                            ind
                        )
                        population[j] = ind
                        gen_scores[j] = after

            # --- SCORE KEEPING --------------------------------+
            gen_avg = sum(gen_scores) / self.popsize
            gen_best = min(gen_scores)

            trial_avg = sum(trial_scores) / self.popsize
            trial_best = min(trial_scores)

            gen_sol = population[gen_scores.index(min(gen_scores))]

            self.logger.info("   > GENERATION AVERAGE: %f " % gen_avg)
            self.logger.info("   > GENERATION BEST: %f " % gen_best)
            self.logger.info(
                "   > TRIAL INFO: {:.2f} {:.2f} ".format(trial_best, trial_avg)
            )

            (
                best_pdb,
                best_SixD_vector,
                best_rmsd,
            ) = self.popul_calculator.cost_func.render_best(i, gen_sol, population)
            best_sol_str = self.popul_calculator.cost_func.get_sol_string(
                best_SixD_vector
            )
            self.logger.info("   > BEST SOL: {} ".format(best_sol_str))
            self.popul_calculator.cost_func.print_information(population)

            # self.popul_calculator.cost_func.pymol_visualization(population)

            name = self.config.out_path + "/evolved.pdb"
            # best_pdb.dump_pdb(name)

            file_object.write("%f \t" % gen_avg)
            file_object.write("%f \t" % gen_best)
            file_object.write("%f \n" % best_rmsd)
            file_object.close()
            gen_end = time.time()
            self.logger.info(f"selection stage in {gen_end - gen_start} seconds")
            file_time.write("%f \n" % (gen_end - gen_start))
            file_time.close()

        return population, best_pdb

    def apply_popul_flexbb(self, population):
        return self.flexbb_swap_operator.apply(population)
