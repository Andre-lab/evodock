#!/usr/bin/env python
# coding: utf-8


# todo: add it as a parameter at the config file
LOW_LIMIT_INIT_DIVERSITY = -1


class SingleProcessPopulCalculator:
    def __init__(self, cost_func, config):
        self.config = config
        self.size = 1
        self.local_search = cost_func.local_search
        self.scfxn = cost_func.local_search.scfxn
        self.cost_func = cost_func

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
            else:
                scored_ind, _, _ = self.local_search.process_individual(ind)
            result_pop.append(scored_ind)
        return result_pop

    def terminate(self):
        # terminate
        pass
