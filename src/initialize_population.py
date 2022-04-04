import random
import glob
import time

from src.genotype_converter import RefineCluspro, generate_genotype
from src.utils import make_trial


class InitializePopulationBuilder:
    def run(self, config, logger, popul_calculator):
        docking_type_option = config.docking_type_option
        if docking_type_option == "Global":
            return InitializePopulation(
                config, logger, popul_calculator
            ).init_population()
        elif docking_type_option == "Local":
            return InitializePopulationLocal(
                config, logger, popul_calculator
            ).init_population()
        elif docking_type_option == "Refine":
            return InitializePopulationRefine(
                config, logger, popul_calculator
            ).init_population()
        elif docking_type_option == "Flexbb":
            return InitializePopulationFlexbb(
                config, logger, popul_calculator
            ).init_population()
        else:
            print(
                "initialize population was not possible due to wrong docking type option"
            )
            exit()


class InitializePopulation:
    # Default initizalize Popualtion for Global docking
    def __init__(self, config, logger, popul_calculator):
        self.config = config
        self.logger = logger
        self.popul_calculator = popul_calculator
        self.popsize = config.popsize
        self.ind_size = 6
        self.bounds = [(-1, 1)] * self.ind_size

    def init_population(self):
        # --- INITIALIZE A POPULATION (step #1) ----------------+
        population_calculator = self.popul_calculator
        popsize = self.popsize
        if popsize is None:
            popsize = self.popsize
        self.logger.info(" init population")
        popul = []
        for i in range(0, popsize):
            indv = []
            for j in range(len(self.bounds)):
                indv.append(random.uniform(self.bounds[j][0], self.bounds[j][1]))

            popul.append(make_trial(i, indv))

        init_population = True
        start = time.time()
        population = population_calculator.run(popul, init_population)
        end = time.time()
        self.logger.info(f" population init in {end - start:.2f} seconds")
        # self.popul_calculator.pymol_visualization(population)
        return population


class InitializePopulationRefine(InitializePopulation):
    def init_population(self):
        population_calculator = self.popul_calculator
        popsize = self.popsize
        refCluspro = RefineCluspro(self.config, self.max_translation)

        if popsize is None:
            popsize = self.popsize
        self.logger.info(" init population")
        popul = []
        for i in range(0, popsize):
            indv = refCluspro.refine_cluspro(i)
            popul.append(make_trial(i, indv))

        init_population = True
        start = time.time()
        population = population_calculator.run(popul, init_population)
        end = time.time()
        self.logger.info(f" population init in {end - start:.2f} seconds")
        # self.popul_calculator.pymol_visualization(population)
        return population


class InitializePopulationLocal(InitializePopulation):
    def init_population(self):
        population_calculator = self.popul_calculator
        popsize = self.popsize
        if popsize is None:
            popsize = self.popsize
        self.logger.info(" init population")
        # default values
        popul = []
        for i in range(0, popsize):
            indv = generate_genotype(self.scfxn.input_pose, self.max_translation)

            popul.append(make_trial(i, indv))

        init_population = True
        start = time.time()
        population = population_calculator.run(popul, init_population)
        end = time.time()
        self.logger.info(f" population init in {end - start:.2f} seconds")
        # self.popul_calculator.pymol_visualization(population)
        return population


class InitializePopulationFlexbb(InitializePopulation):
    def init_population(self):
        population_calculator = self.popul_calculator
        popsize = self.popsize
        if popsize is None:
            popsize = self.popsize
        self.logger.info(" init population")
        popul = []
        # length of flexbb library
        len_receptors = len(glob.glob(self.config.path_receptors)) - 1
        len_ligands = len(glob.glob(self.config.path_ligands)) - 1

        for i in range(0, popsize):
            indv = []
            for j in range(len(self.bounds)):
                indv.append(random.uniform(self.bounds[j][0], self.bounds[j][1]))

            idx_receptor = random.randint(0, len_receptors)
            idx_ligand = random.randint(0, len_ligands)
            popul.append(make_trial(i, indv, idx_ligand, idx_receptor))

        init_population = True
        start = time.time()
        population = population_calculator.run(popul, init_population)
        end = time.time()
        self.logger.info(f" population init in {end - start:.2f} seconds")
        # self.popul_calculator.pymol_visualization(population)
        return population
