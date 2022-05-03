import os

# from pyrosetta.rosetta.protocols.moves import PyMOLMover

# from src.utils import IP_ADDRESS

LOW_LIMIT_INIT_DIVERSITY = 2


class ScorePopulation:
    # def __init__(self, scfxn, jobid, local_search_opt="None", config=None):
    def __init__(self, config, scfxn):
        self.config = config
        self.name = "ScorePopulation"
        self.out_path = config.out_path
        self.scfxn = scfxn
        self.log_best = self.out_path + "/best_individual.csv"
        self.log_interface = self.out_path + "/interface.csv"
        self.log_popul = self.out_path + "/popul.csv"
        self.log_trials = self.out_path + "/trials.csv"
        self.local_search = scfxn.local_search
        self.print_header_logfiles()

    def print_header_logfiles(self):
        with open(self.log_best, "w") as file_object:
            file_object.write("g1,g2,g3,g4,g5,g6\n")
        with open(self.log_trials, "w") as file_object:
            vals = ",".join([f"sc_{t},rmsd_{t}" for t in range(0, self.config.popsize)])
            file_object.write(f"{vals}\n")
        with open(self.log_popul, "w") as file_object:
            vals = ",".join(
                [
                    f"sc_{t},rmsd_{t},Isc_{t},Irmsd_{t}"
                    for t in range(0, self.config.popsize)
                ]
            )
            file_object.write(f"{vals}\n")

    def get_sol_string(self, sol):
        return self.scfxn.get_sol_string(sol)

    def render_best(self, gen, best_solution, population):
        DoFs_vector = self.scfxn.convert_genotype_to_positions(best_solution.genotype)
        rmsd = best_solution.rmsd
        with open(self.log_best, "a") as file_object:
            # file_object.write("gen:\t{}\t".format(gen))
            vector_str = ",".join(["{}".format(i) for i in DoFs_vector])
            file_object.write("{}\n".format(vector_str))

        # best_pdb = self.scfxn.apply_genotype_to_pose(best_solution.genotype)
        best_pdb = self.local_search.best_pose
        return best_pdb, DoFs_vector, rmsd

    def size(self):
        return self.scfxn.size()

    def print_popul_info(self, popul, destiny, trial_popul=False):
        if trial_popul is False:
            popul_str = [
                f"{ind.score:.2f},{ind.rmsd:.2f},{ind.i_sc:.2f},{ind.irms:.2f}"
                for ind in popul
            ]
        else:
            popul_str = [f"{ind.score:.2f},{ind.rmsd:.2f}" for ind in popul]

        with open(destiny, "a") as file_object:
            file_object.write(f"{','.join(popul_str)}\n")

    def print_information(self, popul, trial_popul=False):
        if trial_popul is False:
            destiny = self.log_popul
        else:
            destiny = self.log_trials
        self.print_popul_info(popul, destiny, trial_popul)

    def dump_pdbs(self, popul, destiny="./"):
        os.makedirs(destiny, exist_ok=True)
        for ind, p in enumerate(popul):
            gen = p.genotype
            tmp_pose = self.scfxn.apply_genotype_to_pose(gen)
            tmp_pose.pdb_info().name("popul_" + str(ind))
            tmp_pose.dump_pdb(destiny + "popul_" + str(ind) + ".pdb")

# # <<<<<<< HEAD
#     def pymol_popul_visualization(self, popul, history=False):
#         pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
#         if history:
#             pymover.keep_history(True)
# # =======
#     def pymol_visualization(self, popul):
#         # pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
# # >>>>>>> main
#         for ind, p in enumerate(popul):
#             gen = p.genotype
#             tmp_pose = self.scfxn.apply_genotype_to_pose(gen)
#             tmp_pose.pdb_info().name("popul_" + str(ind))
#             # pymover.apply(tmp_pose)

    def run(self, popul, init_population=False):
        convert_pop = popul
        result_pop = []
        for ind in convert_pop:
            if init_population:
                (
                    scored_ind,
                    before,
                    after,
                ) = self.local_search.process_individual(ind, True)
                if (
                    self.config.docking_type_option == "Refine"
                    and LOW_LIMIT_INIT_DIVERSITY > 0
                ):
                    while scored_ind.rmsd < LOW_LIMIT_INIT_DIVERSITY:
                        (scored_ind, before, after) = self.randomize_ind(scored_ind)
            else:
                scored_ind, _, _ = self.local_search.process_individual(ind)
            result_pop.append(scored_ind)
        return result_pop
