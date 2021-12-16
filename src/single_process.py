#!/usr/bin/env python
# coding: utf-8



from src.individual import Individual
from src.genotype_converter import generate_genotype
# todo: add it as a parameter at the config file
LOW_LIMIT_INIT_DIVERSITY = 2


class SingleProcessPopulCalculator:
    def __init__(self, cost_func, config):
        self.config = config
        self.size = 1
        self.local_search = cost_func.local_search
        self.scfxn = cost_func.local_search.scfxn
        self.cost_func = cost_func
        self.max_translation = config.get_max_translation()

    def randomize_ind(self, scored_ind):
        new_genotype= generate_genotype(
            self.scfxn.input_pose, self.max_translation
        )
        ind = Individual(new_genotype, 1, 1, 0, 1000)
        (
            scored_ind,
            before,
            after,
        ) = self.local_search.process_individual(ind, True)
        return (scored_ind, before, after) 

    def run(self, popul, init_population=False):
        convert_pop = popul
        result_pop = []
        for ind in convert_pop:
            if init_population:
                (
                    scored_ind,
                    before,
                    after,
                ) = self.cost_func.local_search.process_individual(ind, True)
                if LOW_LIMIT_INIT_DIVERSITY > 0:
                    while scored_ind.rmsd < LOW_LIMIT_INIT_DIVERSITY:
                        (scored_ind, before, after) = self.randomize_ind(scored_ind)
            else:
                scored_ind, _, _ = self.local_search.process_individual(ind)
            result_pop.append(scored_ind)
        return result_pop

    def terminate(self):
        # terminate
        pass
