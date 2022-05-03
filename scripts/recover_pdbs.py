""" Get PDB from position and calculate DockQ?"""
#!/usr/bin/env python
# coding: utf-8

import sys
import time
import glob
import configparser
import ast
import os

from pyrosetta import Pose, init
from src.config_reader import EvodockConfig
from src.scfxn_fullatom import FAFitnessFunction
from src.local_search import LocalSearchPopulation
from src.options import init_global_docking
from src.utils import get_pose_from_file

# [inputs]
# pose_input=/easy_dock/1oph_AB.prepack.pdb

# [outputs]
# output_file=./results_LargeDefaultDE/28102020_1oph_default_F09CR03/1oph/4/evolution_f09_cr03.log

# [position]
# rot_max_magnitude=180
# trans_max_magnitude=70

# [DE]
# scheme=BEST
# popsize=100
# mutate=0.9
# recombination=0.3
# maxiter=100
# local_search=mcm_rosetta


data = {
    "1oph": "/home/daniel/projects/Benchmark_Seeds/results_LargeDefaultDE/28102020_1oph_default_F09CR03"
}


MAIN_PATH = "/home/daniel/projects/Benchmark_Seeds/"

data = {
    "1ml0": MAIN_PATH
    + "results_NewConfigLargeDefault/04112020_first4proteins_F09CR03/",
    "2hle": MAIN_PATH
    + "results_NewConfigLargeDefault/04112020_first4proteins_F09CR03/",
    "2hrk": MAIN_PATH + "/results_Large_DefaultDE/20102020_DefaultDE_F09CR03/",
    "1b6c": MAIN_PATH + "/results_Large_DefaultDE/20102020_DefaultDE_More_F09CR03/",
    "1ktz": MAIN_PATH + "/results_Large_DefaultDE/20102020_DefaultDE_More_F09CR03/",
    "1qa9": MAIN_PATH + "/results_Large_DefaultDE/20102020_DefaultDE_More_F09CR03/",
    "1kxp": MAIN_PATH + "/results_LargeDefaultDE/27102020_large_more_F09CR03/",
    "1ppe": MAIN_PATH
    + "results_NewConfigLargeDefault/04112020_first4proteins_F09CR03/",
    "1oph": MAIN_PATH
    + "results_NewConfigLargeDefault/04112020_first4proteins_F09CR03/",
    "1gcq": MAIN_PATH + "./results_LargeDefault/02112020_1cgq_default_F09CR03/",
}


class MockConfig:
    def __init__(self, prot, config):
        evodock_folder = "/home/daniel/projects/evodock/"
        self.docking_type_option = "Global"
        self.bb_strategy = "None"
        self.relax_prob = -1
        self.pose_input = evodock_folder + config["inputs"].get("pose_input")
        self.native_input = self.pose_input
        self.path_ligands = ""
        self.path_receptors = ""
        self.relax_prob = -1
        self.scheme = config["DE"].get("scheme")
        self.popsize = config["DE"].getint("popsize")
        self.mutate = config["DE"].getfloat("mutate")
        self.recombination = config["DE"].getfloat("recombination")
        self.maxiter = config["DE"].getint("maxiter")
        self.local_search_option = "None"
        self.jobid = f"./fake_jobid_{prot}"
        self.output_file = self.jobid + "/evolution.log"

    def get_max_translation(self):
        return 70


def get_config(prot, benchmark_folder):
    conf_path = benchmark_folder + "/configs/*" + prot + "*.ini"
    ini_file = glob.glob(conf_path)[0]
    config = configparser.ConfigParser()
    config.read(ini_file, encoding="utf-8-sig")
    conf = MockConfig(prot, config)
    return conf


def main():

    # prot = "1oph"
    prot = sys.argv[-1]
    benchmark_folder = data[prot]

    config = get_config(prot, benchmark_folder)

    init(extra_options=init_global_docking(config.pose_input))
    config.jobid = f"./fake_jobid_{prot}"
    config.output_file = config.jobid + "/evolution.log"

    native_input = config.pose_input

    input_pose = get_pose_from_file(config.pose_input)
    native_pose = get_pose_from_file(native_input)

    scfxn = FAFitnessFunction(input_pose, native_pose, config)
    # local_search_option = config.local_search_option
    # local_search = LocalSearchPopulation(scfxn, local_search_option, config)

    os.makedirs(f"fast_extract_{prot}", exist_ok=True)
    best_files = glob.glob(benchmark_folder + "/" + prot + "/*/best_f09_cr03.log")
    for idx, bfile in enumerate(best_files):
        with open(bfile, "r") as f:
            for line in f:
                pass
            last_line = line
        positions = list(ast.literal_eval(last_line.strip()))
        pose = scfxn.get_solution_from_positions(positions)
        pose.dump_pdb(f"fast_extract_{prot}/{prot}_{idx}.pdb")


if __name__ == "__main__":
    main()
