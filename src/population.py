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
        self.log_all_geno = self.out_path + "/all_individuals.csv"
        self.log_interface = self.out_path + "/interface.csv"
        self.log_popul = self.out_path + "/popul.csv"
        self.log_trials = self.out_path + "/trials.csv"
        self.log_ensembles = self.out_path + "/ensemble.csv"
        self.log_flipfix = self.out_path + "/flip_fix.csv"
        self.local_search = scfxn.local_search
        self.print_header_logfiles()

    def print_header_logfiles(self):
        with open(self.log_best, "w") as file_object:
            file_object.write("gen,g1,g2,g3,g4,g5,g6\n")
        with open(self.log_all_geno, "w") as file_object:
            if self.config.syminfo:
                genes = self.config.syminfo.genotype_size
            else:
                genes = 6
            vals = f"gen," + ",".join([",".join([f"g{gen}_{ind}" for gen in range(1, genes + 1)])  for ind in range(0, self.config.popsize)])
            file_object.write(f"{vals}\n")
        with open(self.log_trials, "w") as file_object:
            vals = "gen," + ",".join([f"sc_{t},rmsd_{t}" for t in range(0, self.config.popsize)])
            file_object.write(f"{vals}\n")
        with open(self.log_popul, "w") as file_object:
            vals = "gen," + ",".join(
                [
                    f"sc_{t},rmsd_{t},Isc_{t},Irmsd_{t}"
                    for t in range(0, self.config.popsize)
                ]
            )
            file_object.write(f"{vals}\n")
        if self.config.flexbb:
            with open(self.log_ensembles, "w") as file_object:
                file_object.write("gen," + ",".join([str(i) for i in range(0, self.config.popsize)]) + "\n")
            if self.config.docking_type_option == "GlobalFromMultimer":
                with open(self.log_flipfix, "w") as f:
                    f.write('#flipped and fixed state @ all generations\n')

    def get_sol_string(self, sol):
        return self.scfxn.get_sol_string(sol)

    def render_best(self, gen, best_solution, population):
        DoFs_vector = self.scfxn.convert_genotype_to_positions(best_solution.genotype)
        rmsd = best_solution.rmsd
        with open(self.log_best, "a") as file_object:
            # file_object.write("gen:\t{}\t".format(gen))
            vector_str = ",".join(["{}".format(i) for i in DoFs_vector])
            file_object.write(str(gen) + "," + "{}\n".format(vector_str))

        # best_pdb = self.scfxn.apply_genotype_to_pose(best_solution.genotype)
        best_pdb = self.local_search.best_pose
        return best_pdb, DoFs_vector, rmsd

    def print_genotype_values(self, population, gen):
        with open(self.log_all_geno, "a") as file_object:
            vals = f"{gen}," + ",".join([",".join(map(str, self.scfxn.convert_genotype_to_positions(ind.genotype))) for ind in population])
            file_object.write(vals + "\n")

    def size(self):
        return self.scfxn.size()

    def print_popul_info(self, popul, destiny, gen, trial_popul=False):
        if trial_popul is False:
            popul_str = [
                f"{ind.score:.2f},{ind.rmsd:.2f},{ind.i_sc:.2f},{ind.irms:.2f}"
                for ind in popul
            ]
        else:
            popul_str = [f"{ind.score:.2f},{ind.rmsd:.2f}" for ind in popul]

        with open(destiny, "a") as file_object:
            file_object.write(f"{gen},{','.join(popul_str)}\n")

    def print_ensembles(self, population, gen):
        # todo: check if sorting is unnecessary?
        if self.config.flexbb:
            with open(self.log_ensembles, "a") as f:
                f.write(str(gen) + "," + ','.join([ind.subunit_name for ind in population]) + "\n")

    def print_flip_fix_info(self, population):
        if self.config.flexbb and self.config.docking_type_option == "GlobalFromMultimer":
            with open(self.log_flipfix, "a") as f:
                f.write("flipped," + ','.join([str(int(ind.flipped)) for ind in population]) + "\n")
                f.write("fixed," + ','.join([str(int(ind.fixed)) for ind in population]) + "\n")
        else:
            print("nope", self.config.flexbb, self.config.docking_type_option)

    def print_information(self, popul, gen, trial_popul=False):
        if trial_popul is False:
            destiny = self.log_popul
        else:
            destiny = self.log_trials
        self.print_popul_info(popul, destiny, gen, trial_popul)

    def dump_pdbs(self, popul, destiny="./"):
        os.makedirs(destiny, exist_ok=True)
        for ind, p in enumerate(popul):
            gen = p.genotype
            tmp_pose = self.scfxn.apply_genotype_to_pose(gen)
            tmp_pose.pdb_info().name("popul_" + str(ind))
            tmp_pose.dump_pdb(destiny + "popul_" + str(ind) + ".pdb")

    def run(self, popul, init_population=False):
        convert_pop = popul
        result_pop = []
        for ind in convert_pop:
            if init_population:
                (
                    scored_ind,
                    before,
                    after,
                ) = self.local_search.process_individual(ind)
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
