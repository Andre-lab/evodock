import os

from pyrosetta.rosetta.protocols.moves import PyMOLMover

from src.local_search import LocalSearchPopulation
from src.utils import IP_ADDRESS


class ScorePopulation:
    def __init__(self, scfxn, jobid, local_search_opt, config):
        self.config = config
        self.name = "ScorePopulation"
        self.jobid = jobid
        self.scfxn = scfxn
        self.local_search_opt = local_search_opt
        self.log_best = jobid.replace("evolution", "best")
        self.log_interface = jobid.replace("evolution", "interface")
        self.log_popul = jobid.replace("evolution", "popul")
        self.log_trials = jobid.replace("evolution", "trials")
        self.local_search = LocalSearchPopulation(scfxn, local_search_opt, config)
        with open(self.log_best, "w") as file_object:
            file_object.write("#{}\n".format(jobid))
        with open(self.log_trials, "w") as file_object:
            file_object.write("#{}\n".format(jobid))
        with open(self.log_popul, "w") as file_object:
            file_object.write("#{}\n".format(jobid))
        with open(self.log_interface, "w") as file_object:
            file_object.write("#{}\n".format(jobid))

    def get_sol_string(self, sol):
        return self.scfxn.get_sol_string(sol)

    def render_best(self, gen, best_solution, population):
        DoFs_vector = self.scfxn.convert_genotype_to_positions(best_solution.genotype)
        rmsd = best_solution.rmsd
        with open(self.log_best, "a") as file_object:
            # file_object.write("gen:\t{}\t".format(gen))
            vector_str = ",".join(["{}".format(i) for i in DoFs_vector])
            file_object.write("{}\n".format(vector_str))

        best_pdb = self.scfxn.apply_genotype_to_pose(best_solution.genotype)
        return best_pdb, DoFs_vector, rmsd

    def size(self):
        return self.scfxn.size()

    def print_popul_info(self, popul, destiny, trial_popul=False):
        popul_dst = [ind.score for ind in popul]
        popul_rmsd = [ind.rmsd for ind in popul]
        popul_interface = [ind.i_sc for ind in popul]
        popul_irmsd = [ind.irms for ind in popul]

        with open(destiny, "a") as file_object:
            popul_dst_str = ",".join(["{:.2f}".format(i) for i in popul_dst])
            popul_rmsd_str = ",".join(["{:.2f}".format(i) for i in popul_rmsd])
            file_object.write("{}\n".format(popul_dst_str))
            file_object.write("{}\n".format(popul_rmsd_str))

        if trial_popul is False:
            with open(self.log_interface, "a") as file_object:
                popul_interface_str = ",".join(
                    ["{:.2f}".format(i) for i in popul_interface]
                )
                popul_irmsd_str = ",".join(["{:.2f}".format(i) for i in popul_irmsd])
                file_object.write("{}\n".format(popul_interface_str))
                file_object.write("{}\n".format(popul_irmsd_str))

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

    def pymol_visualization(self, popul):
        pymover = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
        for ind, p in enumerate(popul):
            gen = p.genotype
            tmp_pose = self.scfxn.apply_genotype_to_pose(gen)
            tmp_pose.pdb_info().name("popul_" + str(ind))
            pymover.apply(tmp_pose)
