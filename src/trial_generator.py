import random
from src.symmetry import individual_is_within_bounds

from src.utils import make_trial
from operator import itemgetter
from collections import namedtuple


from src.mutation_strategies import (
    MutationStrategyBuilder,
    StrategyRandom,
    StrategyPBest,
    StrategyBest,
    StrategyCurrent,
    StrategyTriangular,
)


class TrialGenerator:
    def __init__(self, config, local_search):
        self.mutation = config.mutate
        self.recombination = config.recombination
        self.mutation_strategy = MutationStrategyBuilder(config).build(size=local_search.scfxn.converter.size)
        self.local_search = local_search
        self.config = config

    def mutate(self, j, population, gen_scores):
        return self.mutation_strategy.create_donor(j, population, gen_scores)

    def recombine(self, j, population, v_donor):
        x_t = population[j].genotype
        v_trial = []
        for k, obj in enumerate(x_t):
            crossover = random.random()
            if crossover <= self.recombination:
                v_trial.append(v_donor[k])
            else:
                v_trial.append(x_t[k])
        return v_trial

    def build(self, j, population, gen_scores):
        v_donor = self.mutate(j, population, gen_scores)
        v_trial = self.recombine(j, population, v_donor)
        ind = make_trial(j, v_trial, 0, 0, 0)
        if self.config.flexbb:
            trial, _, _ = self.local_search.process_individual(ind)
        else:
            trial = self.local_search.process_individual(ind)
        return trial

class FlexbbTrialGenerator(TrialGenerator):
    def __init__(self, config, local_search):
        self.mutation = config.mutate
        self.recombination = config.recombination
        self.mutation_strategy = MutationStrategyBuilder(config).build(size=local_search.scfxn.converter.size)
        self.local_search = local_search

    def build(self, j, population, gen_scores):
        individual_is_within_bounds(self.local_search.config, self.local_search.scfxn, population[j])
        v_donor = self.mutate(j, population, gen_scores)
        v_trial = self.recombine(j, population, v_donor)
        idx_receptor, idx_ligand, idx_subunit = population[j].idx_receptor, population[j].idx_ligand, population[j].idx_subunit
        receptor_name, ligand_name, subunit_name = population[j].receptor_name, population[j].ligand_name, population[j].subunit_name
        flipped, fixed = population[j].flipped, population[j].fixed
        ind = make_trial(j, v_trial, idx_ligand, idx_receptor, idx_subunit,  receptor_name, ligand_name, subunit_name, flipped, fixed)
        individual_is_within_bounds(self.local_search.config, self.local_search.scfxn, ind)
        trial = self.local_search.process_individual(ind)
        return trial

class pBestGenerator(TrialGenerator):
    def __init__(self, config, rgen):
        self.mutation = config.mutate
        self.recombination = config.recombination
        self.mutation_strategy = MutationStrategyBuilder(config).build(None)
        self.rgen = rgen

    def mutate(self, j, population, gen_scores):
        if random.random() <= self.rgen:
            v_donor = StrategyRandom(self.config).create_donor(
                j, population, gen_scores
            )
        else:
            v_donor = StrategyPBest(self.config).create_donor(j, population, gen_scores)
        return v_donor


class CodeGenerator(TrialGenerator):
    def __init__(self, config, local_search):
        self.mutation = config.mutate
        self.recombination = config.recombination
        self.mutation_strategy = MutationStrategyBuilder(config).build(None)
        self.local_search = local_search

    def random_config(self):
        Config = namedtuple("Config", "mutate, recombination")
        comb_1 = [0.1, 0.5]
        comb_2 = [0.1, 0.3]
        comb_3 = [0.5, 0.1]
        combinations = [comb_1, comb_2, comb_3]
        select_comb = combinations[random.randrange(0, len(combinations))]
        return Config(select_comb[0], select_comb[1])

    def build(self, j, population, gen_scores):
        rand_comb_1 = self.random_config()
        v_donor_1 = StrategyRandom(rand_comb_1).create_donor(j, population, gen_scores)
        rand_comb_2 = self.random_config()
        v_donor_2 = StrategyBest(rand_comb_2).create_donor(j, population, gen_scores)
        rand_comb_3 = self.random_config()
        v_donor_3 = StrategyCurrent(rand_comb_3).create_donor(j, population, gen_scores)
        v_trial_1 = self.recombine(j, population, v_donor_1)
        v_trial_2 = self.recombine(j, population, v_donor_2)
        v_trial_3 = v_donor_3
        ind_1 = make_trial(j, v_trial_1, 0, 0)
        ind_2 = make_trial(j, v_trial_2, 0, 0)
        ind_3 = make_trial(j, v_trial_3, 0, 0)
        trial_1, _, score1 = self.local_search.process_individual(ind_1)
        trial_2, _, score2 = self.local_search.process_individual(ind_2)
        trial_3, _, score3 = self.local_search.process_individual(ind_3)
        data = [(trial_1, score1), (trial_2, score2), (trial_3, score3)]
        trial = min(data, key=itemgetter(1))[0]
        return trial


class TriangularGenerator(TrialGenerator):
    def __init__(self, config, rgen, local_search):
        self.config = config
        self.mutation = config.mutate
        self.recombination = config.recombination
        self.mutation_strategy = MutationStrategyBuilder(config).build(None)
        self.rgen = rgen
        self.local_search = local_search
        self.recombination = 0.8 + (0.1 - 0.8) * pow(
            1 - self.rgen / self.config.maxiter, 4
        )
        # self.recombination = 0.1 + ((0.5 - 0.1) / self.config.maxiter) * self.rgen
        # self.recombination = 1.0

    def mutate(self, j, population, gen_scores):
        if random.random() >= pow(1 - (self.rgen / 100), 2):
            # if random.random() <= 0.5:
            v_donor = StrategyTriangular(self.config).create_donor(
                j, population, gen_scores
            )
        else:
            v_donor = StrategyRandom(self.config).create_donor(
                j, population, gen_scores
            )

        return v_donor

    def build(self, j, population, gen_scores):
        v_donor = self.mutate(j, population, gen_scores)
        v_trial = self.recombine(j, population, v_donor)
        ind = make_trial(j, v_trial, 0, 0)
        trial, before, after = self.local_search.process_individual(ind)
        return trial
