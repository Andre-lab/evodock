"""Script to create experiments batch"""
#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys

import pandas as pd
import numpy as np


""" DOCUMENTATION """

parser = argparse.ArgumentParser(
    """
This script creates a folder herarchy for the bencharmk dataset
    
input: 
     (required)
     - csv: dataset with pdbid and path of each protein, in comma separated .csv format
     - exp: name of the experiment
     (optional)
     - pdbid: desired pdbid to be created, by default = all

output:
     Folders under the directory with name {exp}.
     The experiments folder contains a directory for each pdbid in dataset (or single pdbid if its indicated with the --pdbid <name>)
     Under each protein folder experiment, a folder for RosettaDOCK and EvoDOCK
     Under each type of algorithm, the command lines and the config files if neccesary
     A table information about experiment folder
     print at screen the command lines to execute


examples:
     create dataset.csv:
     pdbid,path
     1b6c,inputs/pdbs/1b6c_AB.prepack.pdb 
     1b6c,inputs/pdbs/1ktz_AB.prepack.pdb

    
     * to view the help:
     python experiments_creator.py -h
     * to run with all complex in the dataset and create the experiment <name>:
     python experiments_creator.py --csv inputs/dataset.csv --exp experiment_Feb11
     * to run with pdbid complex in the dataset and create the experiment <name>:
     python experiments_creator.py --csv inputs/dataset.csv --exp experiment_Feb11 --pdbid 1k6z


requirements:
     The pdbs should be clean (clean_pdb.py) and prepacked. 
     NOTE: The dataset.csv should contain the paths to the '.prepack.pdb' of each complex!
           The prepack_cmd.sh is an example of how to run it, but it should have been used before.
"""
)

""" ARGUMETNS INPUTS """

parser.add_argument("--csv", help="dataset with pdbs in .csv format", required=True)
parser.add_argument("--exp", help="name of the experiment", required=True)
parser.add_argument("--pdbid", help="pdbid from dataset to be created", default="all")
args = parser.parse_args()


""" AUXILIARY FUNCTIONS """


def make_rosetta_experiment(protein_exp, row):
    rosetta_exp = f"{protein_exp}/RosettaDOCK"
    os.makedirs(rosetta_exp, exist_ok=True)
    # with open(f"{rosetta_exp}/prepack_cmd.sh", "w") as f:
    #     f.write("dockingprotocol.linuxgccrelease /\n")
    #     f.write(f"−s {row.path.replace('.prepack','')} /\n")
    #     f.write("−partners A_B /\n")
    #     f.write("−dockppk / \n")
    #     f.write("−ex1 −ex2aro / \n")
    #     f.write("−extrachicutoff 0 /\n")
    #     f.write(f"−unboundrot {row.path.replace('.prepack','')}\n")

    with open(f"{rosetta_exp}/docking_cmd.sh", "w") as f:
        f.write(
            "/home/ingemar/git/Rosetta_shape/main/source/bin/dockingprotocol.linuxgccrelease /\n"
        )
        f.write(f"−s {row.path} /\n")
        f.write("−partners A_B /\n")
        f.write("−nstruct 1/ \n")
        f.write("−randomize1 / \n")
        f.write("−randomize2 / \n")
        f.write("−spin / \n")
        f.write("−ex1 −ex2aro / \n")
        f.write("−extrachicutoff 0 /\n")
        f.write(f"−unboundrot {row.path.replace('.unbound', '')}\n")


def make_evodock_experiment(protein_exp, row):
    evodock_exp = f"{protein_exp}/EvoDOCK"
    os.makedirs(evodock_exp, exist_ok=True)

    execs = []
    for i in range(100):
        config_name = f"{evodock_exp}/config_{row.pdbid}_{i}.ini"
        with open(config_name, "w") as f:
            f.write("[Docking]\n")
            f.write("type=Global\n")
            f.write("bb_strategy=None\n")
            f.write("[inputs]\n")
            f.write(f"pose_input={row.path}\n")
            f.write(f"native_input={row.path.replace('.unbound', '')}\n")
            f.write("path_ligands=''\n")
            f.write("path_receptors=''\n")
            f.write("[outputs]\n")
            f.write(f"output_file={evodock_exp}_{i}/quick_evolution_sample.log\n")
            f.write("[DE]\n")
            f.write("scheme=BEST\n")
            f.write("popsize=100\n")
            f.write(f"mutate={mutacion}\n")
            f.write(f"recombination={crossover}\n")
            f.write("maxiter=100\n")
            f.write("local_search=mcm_rosetta\n")

        exec = f"python evodock.py {config_name} &\n"
        execs.append(exec)
        # with open(f"{evodock_exp}/cmd_{i}.sh", "w") as f:
        #    f.write(f"python evodock.py {config_name}\n")

    return execs


def make_table_experiment(exp, df):
    with open(f"{exp}/table.md", "w") as f:
        f.write("|pdbid|path|\n")
        for idx, row in df.iterrows():
            f.write(f"|{row.pdbid}|{row.path}|\n")

    print(f"2. create informative table at {exp}/table.md")


def read_dataset(args):
    dataset = pd.read_csv(args.csv)

    if args.pdbid != "all":
        dataset = dataset[dataset.pdbid == args.pdbid]
        if len(dataset) == 0:
            print(f"the pdbid {args.pdbid} cannot be found in the dataset")
            exit()

    complex_list = {",".join(dataset.pdbid.to_list())}
    print(f"0. read dataset with {len(dataset)} complexes: {complex_list}")
    return dataset


def print_final_information_and_commandlines(dataset):
    print(f"3. folders herarchy creaded for {','.join(dataset.pdbid.to_list())}")

    print("\n")
    print("===== command lines ========")
    for idx, row in dataset.iterrows():
        protein_exp = f"{args.exp}/{row.pdbid}"
        evodock_exp = f"sh {protein_exp}/EvoDOCK/cmd.sh"
        rosetta_exp = f"sh {protein_exp}/EvoDOCK/docking_cmd.sh"
        # print(evodock_exp)
        # print(rosetta_exp)


def make_directories_for_experiment(args, dataset):
    comp_time = dict(
        zip(dataset.pdbid.unique(), ["11:00:00"] * len(dataset.pdbid.unique()))
    )
    comp_time["1WQ1"] = "7:00:00"
    comp_time["2PCC"] = "7:00:00"
    comp_time["1WEJ"] = "20:00:00"
    os.makedirs(args.exp, exist_ok=True)
    for idx, row in dataset.iterrows():
        # make_rosetta_experiment(protein_exp, row)
        for f in ["0.15", "0.20", "0.25", "0.35", "0.40", "0.45"]:
            for cr in ["0.3", "0.7", "0.9"]:
                # for f in ["0.3", "0.7", "0.9"]:
                #     for cr in ["0.3", "0.7", "0.9"]:
                if np.random.randint(100) < 100:
                    protein_exp = f"{args.exp}/{row.pdbid}/{f}/{cr}"
                    os.makedirs(protein_exp, exist_ok=True)
                    execs = make_evodock_experiment(protein_exp, row, f, cr)
                    with open(f"{protein_exp}/batch_{row.pdbid}_{f}_{cr}.sh", "w") as f:
                        f.write("#!/bin/bash\n")
                        f.write(f"#SBATCH -t {comp_time[row.pdbid]}\n")
                        f.write(f"#SBATCH -o output_{row.pdbid}.log\n")
                        f.write(f"#SBATCH -e error_{row.pdbid}.log\n")
                        f.write("#SBATCH -A SNIC2021-5-351\n")
                        f.write("#SBATCH -n 100\n")
                        f.write("#SBATCH --mem-per-cpu=9000MB\n")
                        f.write(
                            "export PYTHONPATH=$PYTHONPATH:/home/d/dvarela/pfs/PyRosetta4.MinSizeRel.python37.linux.release-265/\n"
                        )
                        f.write(
                            "export PYTHONPATH=$PYTHONPATH:/home/d/dvarela/store/evodock\n"
                        )
                        f.write("source ~/store/venv/bin/activate\n")
                        f.write("module load GCCcore/8.2.0 Python/3.7.2\n")
                        for exec in execs:
                            f.write(exec)
                        f.write("wait\n")


def initial_checks(args, dataset):
    print(f"1. initial checks with dataset {args.csv}")
    for idx, row in dataset.iterrows():
        if os.path.isfile(row.path):
            print(f"cannot find file {row.path}")


def main():
    """MAIN FLOW OF THE SCRIPT"""
    dataset = read_dataset(args)
    initial_checks(args, dataset)
    make_directories_for_experiment(args, dataset)
    make_table_experiment(args.exp, dataset)
    print_final_information_and_commandlines(dataset)


if __name__ == "__main__":
    main()
