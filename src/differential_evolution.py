import glob
import logging
import random
import time

from src.genotype_converter import generate_genotype
from src.individual import Individual
from src.mpi_utils import IndividualMPI
from src.selection import GreedySelection
from src.single_process import SingleProcessPopulCalculator


def ensure_bounds(vec, bounds):
    vec_new = []
    # cycle through each variable in vector
    for i, k in enumerate(vec):

        # variable exceedes the minimum boundary
        if vec[i] < bounds[i][0]:
            vec_new.append(bounds[i][0])

        # variable exceedes the maximum boundary
        if vec[i] > bounds[i][1]:
            vec_new.append(bounds[i][1])

        # the variable is fine
        if bounds[i][0] <= vec[i] <= bounds[i][1]:
            vec_new.append(vec[i])

    return vec_new


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
        self.job_id = config.jobid
        self.ind_size = popul_calculator.cost_func.size()
        self.bounds = [(-1, 1)] * self.ind_size
        self.file_time_name = self.job_id.replace("evolution", "time")
        self.max_translation = config.get_max_translation()
        self.init_file()

    def init_population(self, popsize=None, docking_type="Global"):
        # --- INITIALIZE A POPULATION (step #1) ----------------+

        # population_calculator = SingleProcessPopulCalculator(
        #     self.popul_calculator.cost_func, self.config
        # )

        population_calculator = self.popul_calculator

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
                indv = generate_genotype(
                    self.popul_calculator.scfxn.input_pose, self.max_translation
                )

            # todo: read file and count *.pdb
            idx_receptor = random.randint(1, 10)
            idx_ligand = random.randint(1, 10)

            popul.append(make_trial(i, indv, idx_ligand, idx_receptor))
            population.append(Individual(indv, idx_ligand, idx_receptor, 0, 1000))

        init_population = True
        population = population_calculator.run(popul, init_population)
        self.popul_calculator.cost_func.pymol_visualization(population)
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
        for i in range(1, self.maxiter + 1):
            self.logger.info(" GENERATION:" + str(i))
            start = time.time()
            file_object = open(self.job_id, "a")
            file_time = open(self.file_time_name, "a")

            file_object.write("GENERATION: \t" + str(i) + "\t")
            gen_scores = [ind.score for ind in population]

            gen_sol = population[gen_scores.index(min(gen_scores))]
            # cycle through each individual in the population
            trials = []
            for j in range(0, self.popsize):

                # --- MUTATION (step #3.A) ---------------------+
                # select 3 random vector index positions [0, self.popsize)
                # not including current vector (j)
                candidates = list(range(0, self.popsize))
                candidates.remove(j)
                random_index = random.sample(candidates, 3)

                if self.scheme == "CURRENT":
                    x_1 = population[j].genotype
                if self.scheme == "RANDOM":
                    x_1 = population[random_index[0]].genotype
                if self.scheme == "BEST":
                    x_1 = population[gen_scores.index(min(gen_scores))].genotype

                x_2 = population[random_index[1]].genotype
                x_3 = population[random_index[2]].genotype
                x_t = population[j].genotype  # target individual

                # subtract x3 from x2, and create a new vector (x_diff)
                x_diff = [x_2_i - x_3_i for x_2_i, x_3_i in zip(x_2, x_3)]

                # multiply x_diff by the mutation factor (F) and add to x_1
                v_donor = [
                    x_1_i + self.mutate * x_diff_i
                    for x_1_i, x_diff_i in zip(x_1, x_diff)
                ]
                v_donor = ensure_bounds(v_donor, self.bounds)

                # --- RECOMBINATION (step #3.B) ----------------+

                v_trial = []
                for k, obj in enumerate(x_t):
                    crossover = random.random()
                    if crossover <= self.recombination:
                        v_trial.append(v_donor[k])
                    else:
                        v_trial.append(x_t[k])
                trials.append(v_trial)

            # --- SELECTION (step #3.C) -------------+
            trials = [
                make_trial(idx, t, gen_sol.idx_ligand, gen_sol.idx_receptor)
                for idx, t in enumerate(trials)
            ]
            trial_inds = self.popul_calculator.run(trials)
            self.popul_calculator.cost_func.print_information(trial_inds, True)

            population, gen_scores, trial_scores = GreedySelection().apply(
                trial_inds, population
            )

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
            self.popul_calculator.cost_func.pymol_visualization(population)

            file_object.write("%f \t" % gen_avg)
            file_object.write("%f \t" % gen_best)
            file_object.write("%f \n" % best_rmsd)
            file_object.close()
            end = time.time()
            file_time.write("%f \n" % (end - start))
            file_time.close()

        return population, best_pdb
