import copy
import time

import numpy as np

from src.differential_evolution import Individual

# from mpi4py import MPI


# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()


class IndividualMPI:
    def __init__(self, idx, ind):
        self.idx = idx
        if type(ind) is list:
            self.genotype = ind
            self.score = 1000
            self.rmsd = 0
            self.i_sc = 0
            self.irms = 0
        else:
            self.genotype = ind.genotype
            self.score = ind.score
            self.rmsd = ind.rmsd
            self.i_sc = ind.i_sc
            self.irms = ind.irms

    def convert_to_np(self):
        a = self.genotype
        a.append(self.score)
        a.append(self.rmsd)
        a.append(self.i_sc)
        a.append(self.irms)
        a.append(self.idx)
        a = np.array(a)
        return a

    def convert_to_ind(self):
        return Individual(self.genotype, self.score, self.rmsd, self.i_sc, self.irms)


def np_to_ind(a):
    a = list(a)
    ind = Individual(a[0:6], a[6], a[7], a[8], a[9])
    idx = int(round(a[-1]))
    return (idx, ind)


class MasterProcess:
    def __init__(self, size, cost_func):
        self.rank = 0
        self.size = size
        self.cost_func = cost_func

    def __make_chunks(self, popul, size):
        lst = []
        for d in range(0, len(popul)):
            lst.append(IndividualMPI(d, popul[d]).convert_to_np())
        size = size
        chunks = [lst[i : i + size] for i in range(0, len(lst), size)]
        return chunks

    def run(self, popul):
        list_popul = popul
        req = []
        ind_per_process = int(len(list_popul) / self.size)
        data = self.__make_chunks(list_popul, ind_per_process)
        ind_size = len(data[0][0])
        recv_data = np.zeros(
            shape=(self.size, ind_per_process, ind_size), dtype=np.float64
        )
        run_score_function = 1
        for d in range(1, self.size):
            comm.send(run_score_function, dest=d)
            comm.send(data[d], dest=d)
            r = comm.irecv(buf=recv_data[d], source=d, tag=d)
            req.append(r)

        convert_pop = [np_to_ind(x) for x in data[0]]
        mpi_pop = []
        for idx, ind in convert_pop:
            scored_ind, _, _ = self.cost_func.local_search.process_individual(
                ind.genotype
            )
            mpi_pop.append(IndividualMPI(idx, scored_ind).convert_to_np())

        for idx, ind in enumerate(mpi_pop):
            recv_data[0][idx] = ind

        all_finish = False
        while all_finish is False:
            check_finish = []
            for r in req:
                s = r.Test()
                check_finish.append(s)
            all_finish = all([s is True for s in check_finish])

        # print("recv data ", recv_data)
        result_pop = popul
        for proc in recv_data:
            for arr in proc:
                idx, ind = np_to_ind(arr)
                result_pop[idx] = ind

        return result_pop

    def terminate(self):
        terminate_function = 3
        for d in range(1, self.size):
            comm.send(terminate_function, dest=d)


class Worker:
    def __init__(self, rank, size, cost_func):
        self.rank = rank
        self.size = size
        self.cost_func = cost_func

    def run(self):
        # functions = {1: "run_score", 2: "run_score_rmsd", 3: "terminate"}
        while True:
            action = comm.recv(source=0)
            if action == 3:
                break
            data = comm.recv(source=0)
            # process data here
            convert_pop = [np_to_ind(x) for x in data]
            mpi_pop = []
            for idx, ind in convert_pop:
                scored_ind, _, _ = self.cost_func.local_search.process_individual(
                    ind.genotype
                )
                mpi_pop.append(IndividualMPI(idx, scored_ind).convert_to_np())

            send_data = np.array(mpi_pop, dtype=np.float64)
            comm.Send(send_data, dest=0, tag=rank)
