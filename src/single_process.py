#!/usr/bin/env python
# coding: utf-8


import random

import numpy as np

from src.local_search import LocalSearchPopulation
from src.mpi_utils import IndividualMPI, np_to_ind
from src.scfxn_fullatom import FAFitnessFunction


class SingleProcessPopulCalculator:
    def __init__(self, cost_func, syminfo: dict = None):
        self.size = 1
        self.cost_func = cost_func
        self.local_search = cost_func.local_search
        self.scfxn = cost_func.local_search.scfxn
        scfxn = FAFitnessFunction(self.scfxn.native_pose, self.scfxn.trans_max_magnitude, syminfo)
        init_popul = LocalSearchPopulation(scfxn, "mcm_rosetta")
        self.init_population = init_popul

    def __make_chunks(self, popul, size):
        lst = []
        for d in range(0, len(popul)):
            lst.append(IndividualMPI(d, popul[d]).convert_to_np())
        size = size
        chunks = [lst[i : i + size] for i in range(0, len(lst), size)]
        return chunks

    def run(self, popul, init_population=False):
        list_popul = popul
        ind_per_process = int(len(list_popul) / self.size)
        data = self.__make_chunks(list_popul, ind_per_process)
        ind_size = len(data[0][0])
        recv_data = np.zeros(
            shape=(self.size, ind_per_process, ind_size), dtype=np.float64
        )
        convert_pop = [np_to_ind(x, genotype_size=self.scfxn.size()) for x in data[0]]
        mpi_pop = []
        for idx, ind in convert_pop:
            if ind.score == 1000:
                continue_score = True
                while continue_score is True:
                    if init_population:
                        (
                            scored_ind,
                            before,
                            after,
                        ) = self.cost_func.local_search.process_individual(
                            ind.genotype, True
                        )
                    else:
                        scored_ind, _, _ = self.local_search.process_individual(
                            ind.genotype
                        )

                    if init_population:
                        if scored_ind.rmsd > 2:
                            continue_score = False
                            mpi_pop.append(
                                IndividualMPI(idx, scored_ind).convert_to_np()
                            )
                        else:
                            indg = []
                            for j in range(self.scfxn.size()):
                                indg.append(random.uniform(-1, 1))
                            ind.genotype = indg
                            ind.score = 1000
                    else:
                        continue_score = False
                        mpi_pop.append(IndividualMPI(idx, scored_ind).convert_to_np())
            else:
                mpi_pop.append(IndividualMPI(idx, ind).convert_to_np())

        for idx, ind in enumerate(mpi_pop):
            recv_data[0][idx] = ind

        result_pop = popul
        for proc in recv_data:
            for arr in proc:
                idx, ind = np_to_ind(arr, genotype_size=self.scfxn.size())
                result_pop[idx] = ind

        return result_pop

    def terminate(self):
        # terminate
        pass
