import random
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
import glob
import time

from src.genotype_converter import RefineCluspro, generate_genotype
from src.utils import make_trial
from src.symmetry import individual_is_within_bounds


class InitializePopulationBuilder:
    def run(self, DE):
        docking_type_option = DE.config.docking_type_option
        if docking_type_option == "Global":
            return InitializePopulation(
                DE.config, DE.logger, DE.popul_calculator, DE.scfxn
            ).init_population()
        elif docking_type_option == "Local":
            return InitializePopulationLocal(
                DE.config, DE.logger, DE.popul_calculator, DE.scfxn
            ).init_population()
        elif docking_type_option == "Refine":
            return InitializePopulationRefine(
                DE.config, DE.logger, DE.popul_calculator, DE.scfxn
            ).init_population()
        elif docking_type_option == "Flexbb":
            return InitializePopulationFlexbb(
                DE.config, DE.logger, DE.popul_calculator, DE.scfxn
            ).init_population()
        else:
            print(
                "initialize population was not possible due to wrong docking type option"
            )
            exit()


class InitializePopulation:
    # Default initizalize Popualtion for Global docking
    def __init__(self, config, logger, popul_calculator, scfxn):
        self.config = config
        self.logger = logger
        self.popul_calculator = popul_calculator
        self.popsize = config.popsize
        if config.syminfo:
            self.ind_size = config.syminfo.genotype_size
        else:
            self.ind_size = 6
        self.bounds = [(-1, 1)] * self.ind_size
        self.scfxn = scfxn

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
        return population


class InitializePopulationRefine(InitializePopulation):
    def init_population(self):
        population_calculator = self.popul_calculator
        popsize = self.popsize
        refCluspro = RefineCluspro(self.config, self.config.get_max_translation())

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
            if is_symmetric(self.scfxn.dock_pose):
                # todo:
                #  generate_genotype just sets the local bounds (8 degrees rotation + 3 Å translation) for the dimeric case.
                #  for symmetry i still choose to use the symbounds set in the config file. I dont have default values as of yet: todo!
                indv = []
                if self.config.syminfo.init_bounds:
                    bounds = self.config.syminfo.init_bounds
                else:
                    bounds = self.bounds
                for j in range(len(bounds)):
                    indv.append(random.uniform(bounds[j][0], bounds[j][1]))
            else:
                indv = generate_genotype(self.scfxn.dock_pose, self.max_translation)
            popul.append(make_trial(i, indv))
        init_population = True
        start = time.time()
        population = population_calculator.run(popul, init_population)
        end = time.time()
        self.logger.info(f" population init in {end - start:.2f} seconds")
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
        if self.config.syminfo:
            len_subunits = len(glob.glob(self.config.path_subunits)) - 1
        else:
            len_receptors = len(glob.glob(self.config.path_receptors)) - 1
            len_ligands = len(glob.glob(self.config.path_ligands)) - 1
            self.logger.info(
                f" flexbb with {len_ligands} ligands and {len_receptors} receptors"
            )
        for i in range(0, popsize):
            indv = []
            for j in range(len(self.bounds)):
                indv.append(random.uniform(self.bounds[j][0], self.bounds[j][1]))
            if self.config.syminfo:
                idx_subunit = random.randint(0, len_subunits)
                popul.append(make_trial(i, indv, subunit=idx_subunit))
            else:
                idx_receptor = random.randint(0, len_receptors)
                idx_ligand = random.randint(0, len_ligands)
                popul.append(make_trial(i, indv, ligand=idx_ligand, receptor=idx_receptor))
        init_population = True
        start = time.time()
        for ind in popul:
            individual_is_within_bounds(self.config, self.scfxn, ind)
        population = population_calculator.run(popul, init_population)
        for ind in population:
            individual_is_within_bounds(self.config, self.scfxn, ind)
        end = time.time()
        self.logger.info(f" population init in {end - start:.2f} seconds")
        # add the current names to each individual
        for ind in population:
            if self.config.syminfo:
                ind.subunit_name = self.flexbb_swap_operator.swap_operator.list_subunits[ind.idx_subunit]
            else:
                raise NotImplementedError
        return population
