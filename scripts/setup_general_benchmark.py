#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sets up the benchmark
@Author: Mads Jeppesen
@Date: 7/14/21
"""

from pathlib import Path
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
import textwrap
from itertools import product
import argparse
from distutils.util import strtobool
import shutil
import pandas as pd
from shapedesign.settings import SYMMETRICAL
import numpy as np
import random
from pyrosetta import pose_from_file, init
from af_to_evodock_ensemble import AF2EvoDOCK

# create config files and change the output file
ROOT = Path("..").resolve()
BENCHMARK = ROOT.joinpath("benchmark")

# def create_nma_script(fin):
#     with open(fin, "w") as fout:
#         fout.write("""<ROSETTASCRIPTS>
#     <SCOREFXNS>
#         <ScoreFunction name="bn15_cart" weights="beta_nov15_cart" />
#     </SCOREFXNS>
#     <RESIDUE_SELECTORS>
#     </RESIDUE_SELECTORS>
#     <TASKOPERATIONS>
#     </TASKOPERATIONS>
#     <FILTERS>
#     </FILTERS>
#     <MOVERS>
#         <NormalModeRelax name="nma" cartesian="true" centroid="false"
# scorefxn="bn15_cart" nmodes="5" mix_modes="true" pertscale="1.0"
# randomselect="false" relaxmode="relax" nsample="1"
# cartesian_minimize="false" />
#     </MOVERS>
#     <APPLY_TO_POSE>
#     </APPLY_TO_POSE>
#     <PROTOCOLS>
#         <Add mover="nma" />
#     </PROTOCOLS>
#     <OUTPUT scorefxn="bn15_cart" />
# </ROSETTASCRIPTS>
# """)

def create_node_file(fin, mailer=0, mimir=0, mogwai=0):
    with open(fin, "w") as fout:
        fout.write(f"""{mailer}/:
{mogwai}/mogwai
{mimir}/mimir
#mads@130.235.135.47
#mads@130.235.135.80
#mads@130.235.135.23""")

def create_parallel_input_file(fin, input_files):
    with open(fin, "w") as fout:
        fout.write("\n".join(map(str, input_files)))

def create_parallel_script(parallel_script_path, command, nodefile, joblog, parallel_input_file):
    assert Path(nodefile).exists()
    with open(parallel_script_path, "w") as fout:
        fout.write(f"""#!/bin/bash
parallel_='parallel --link --delay .2 --joblog {joblog} --resume-failed --ungroup --sshloginfile {nodefile}'
command_='{command} :::: {parallel_input_file}'
echo $parallel_ $command_
$parallel_ $command_""")

def create_from_flexbb_template(native_input, native_symmetric_file, native_symdef_file, input_symdef_path, output_path, config_name,
                                subunits, bounds, popsize=100, mutate=0.1, recombination=0.7, maxiter=100, slide=True, low_memory=False,
                                normalize_score=False, scheme="RANDOM", selection="interface", init=None, allow_flip=False,
                                init_input_fix_percent=0):

    with open(config_name, "w") as f:
        f.write(textwrap.dedent(f"""[Docking]
type=GlobalFromMultimer

[Inputs]
subunits={subunits}/
symdef_file={input_symdef_path}

[Native]
crystallic_input={native_input}
symmetric_input={native_symmetric_file}
symdef_file={native_symdef_file}

[Outputs]
output_path={output_path}
output_pdb=True

[Flexbb]
swap_prob=0.3
normalize_score=false
{'low_memory_mode=true' if low_memory else''}


[Bounds]
bounds={bounds}
{'init=' + init if init else''}
initialize_rigid_body_dofs=true
init_input_fix_percent={init_input_fix_percent}
allow_flip={allow_flip}

[Pymol]
on=false
history=false
show_local_search=false

[DE]
scheme={scheme}
popsize={popsize}
mutate={mutate}
recombination={recombination}
maxiter={maxiter}
local_search=symshapedock
slide={slide}
selection={selection}
max_slide_attempts=100

# Initialize Rosetta with these flags
[RosettaOptions]
initialize_rigid_body_dofs=true
"""))

def create_evodock_submission_script_kebnekaise(file_name, ini_files, time="72:00:00", exclusive=True, old_behavior=False):
    with open(file_name, "w") as f:
        f.write(textwrap.dedent(f"""        #!/bin/sh
        #SBATCH -t {time} 
        #SBATCH -o output.log
        #SBATCH -e error.log
        #SBATCH -A SNIC2021-5-351 
        #SBATCH -n {len(ini_files)}
        {"#SBATCH --exclusive" if exclusive else ""}
        #SBATCH --mail-user=mads.jeppesen@biochemistry.lu.se
        #SBATCH --mail-type=END
        module purge > /dev/null 2>&1
        source ~/.bashrc
        conda activate evodock
        """))
        # shorter way that could be sexier
        # f.write(textwrap.dedent(f"""
        #         for ((i=0; i<{len(ini_files)}; i++))
        #         do
        #             srun -Q --exclusive --overlap -n 1 -N 1 evodock.py {ini_files[0]} &> log.txt &
        #             sleep 1
        #         done
        #         """))
        # f.write("wait")
        for ini_file in ini_files:
            f.write(f"srun -N 1 -n 1 --exclusive --overlap python ~/projects/evodock{'2' if old_behavior else ''}/evodock.py {ini_file} &> log.txt &\n")
        f.write("wait")

def create_evodock_submission_script_lunarc(file_name, ini_files, jobname, time="48:00:00", exclusive=True):
    with open(file_name, "w") as f:
        f.write(textwrap.dedent(f"""        #!/bin/sh
        #SBATCH -J {jobname}
        #SBATCH -t {time} 
        #SBATCH -o output.log
        #SBATCH -e error.log
        #SBATCH -A LU2022-2-9
        #SBATCH -n {len(ini_files)}
        {"#SBATCH --exclusive" if exclusive else ""}
        #SBATCH --mail-user=mads.jeppesen@biochemistry.lu.se
        #SBATCH --mail-type=END
        module purge > /dev/null 2>&1
        source ~/.bashrc
        conda activate evodock
        """))
        for n, ini_file in enumerate(ini_files): #THEY ARE SORTED DURING INITIALIZATION
            # get the number
            f.write(f"srun -N 1 -n 1 --exclusive --overlap python ~/projects/evodock/evodock.py configs/{ini_file} &> logs/log{n}.txt &\n")
        f.write("wait")

binpath = "/opt/Rosetta/main/source/bin"
relaxbin = f"{binpath}/relax.linuxgccrelease"
nmabin = f"{binpath}/rosetta_scripts.default.linuxgccrelease "
backrupbin = f"{binpath}/backrub.default.linuxgccrelease "
def create_flexbb_creater_script(preperation_path, subunits_path, native_file_symmetric_cut, runs, ensemble_size):
    # create paths
    scripts_path = preperation_path.joinpath("scripts")
    scripts_path.mkdir(exist_ok=True)
    scores_path = preperation_path.joinpath("scores")
    scores_path.mkdir(exist_ok=True)
    tmp_path = preperation_path.joinpath("tmp")
    tmp_path.mkdir(exist_ok=True)
    # create nma script
    nmascript = scripts_path.joinpath("nma.xml")
    if not nmascript.exists():
        create_nma_script(nmascript)
    # generate the run command
    commands = []
    end_line = f" -in:file:s {{}} -out:suffix {{}}  -overwrite"
    it = iter(range(ensemble_size))
    command_type, suffixes = [], []
    # relax_w_constrains
    for model_name in runs["relax_w_constrains"]:
        suffix = f"__{next(it)}__relaxwconst_"
        suffixes.append(suffix)
        commands.append(f"{relaxbin}  -relax:thorough  -relax:constrain_relax_to_start_coords" + end_line.format(model_name, suffix))
        command_type.append("relaxwconst")
    for model_name in runs["relax"]:
        suffix = f"__{next(it)}__relax_"
        suffixes.append(suffix)
        commands.append(f"{relaxbin}  -relax:thorough " + end_line.format(model_name, suffix))
        command_type.append("relax")
    for model_name in runs["nma"]:
        suffix = f"__{next(it)}__nma_"
        suffixes.append(suffix)
        commands.append(f"{nmabin} -parser:protocol {nmascript} " + end_line.format(model_name, suffix))
        command_type.append("nma")
    for model_name in runs["backrub"]:
        suffix = f"__{next(it)}__backrub_"
        suffixes.append(suffix)
        commands.append(f"{backrupbin} -backrub:ntrials 20000 -backrub:mc_kt 0.6 " + end_line.format(model_name, suffix))
        command_type.append("backrub")
    model_names = runs["relax_w_constrains"] + runs["relax"] + runs["nma"] + runs["backrub"]

    # make script
    scripts = []
    for n, (command, command_type, model_name, suffix) in enumerate(zip(commands, command_type, model_names, suffixes), 1):
        script = scripts_path.joinpath(f"flexbb_script_{command_type}_{n}.py")
        scripts.append(script)
        with open(script, "w") as f:
            f.write(f"""import subprocess
import shutil
import pandas as pd
from pathlib import Path
from pyrosetta import pose_from_file, init
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreFunctionFactory
from shapedesign.src.utilities.alignment import tmalign
import os
init()
# make a tmp directory
tmp_path = Path('{tmp_path.joinpath("tmp" + suffix)}')
tmp_path.mkdir(exist_ok=True) 
os.chdir(tmp_path)
# call the bbsearch script
subprocess.Popen(f'{command} -out:path:pdb {{tmp_path}}"'.split()).wait()
result_name = tmp_path.joinpath('{Path(model_name).stem}' + f'{suffix}_0001.pdb')
if {'True' if not "relax" in command_type else 'False'}: # True == means this is a NMA / backrub run
    print("Found {command_type}:", "the structure will be relaxed.")
    subprocess.Popen(f'{relaxbin} -in:file:s {{result_name}} -relax:thorough -relax:constrain_relax_to_start_coords'.split()).wait()
    subprocess.Popen(f"rm {{result_name}}".split()).wait()
    result_name = tmp_path.joinpath('{Path(model_name).stem}' + f'{suffix}_0001_0001.pdb')
# If the structure is either nma or backrub we have to relax them
################
# ALIGN SCRIPT #
################
# read in 
rot_pose = pose_from_file(str(result_name))
fix_pose = pose_from_file('{native_file_symmetric_cut}')
assert rot_pose.size() == fix_pose.size(), "fix and rot pose are not of equal length."
# align
tmscore = tmalign(rot_pose, fix_pose)
rmsd = CA_rmsd(fix_pose, rot_pose)
# align to the native and check the rmsd
final_name = '{subunits_path}/{Path(model_name).stem}_{command_type}_{n}.pdb'
rot_pose.dump_pdb(final_name)
sfxn = ScoreFunctionFactory.create_score_function("ref2015")
score = sfxn.score(rot_pose)
# output the total score and the rmsd
pd.DataFrame({{'score': [score], 'rmsd':[rmsd], 'tmscore':[tmscore]}} ).to_csv(f'{scores_path}/{{Path(final_name).stem}}' + '.sc', index=False)
subprocess.Popen(f"rm -r {{tmp_path}}".split()).wait()
subprocess.Popen("rm score__*".split()).wait()
subprocess.Popen("rm *last.pdb *low.pdb".split()).wait()

""")
    return scripts


def make_cut_and_copy_script(native_file_symmetric, pdbid, meta_file):
    out = native_file_symmetric.parent.parent.joinpath("preperation/cut_script.py")
    with open(out, "w") as f:
        f.write(f"""from pathlib import Path
from pyrosetta import pose_from_file, init
import pandas as pd
native_relaxed_path = Path('{native_file_symmetric}')
pdbid = '{pdbid}'
meta_file = '{meta_file}'
# read in
init()
full_name = native_relaxed_path.parent.joinpath(f"{{pdbid}}.cif")
pose = pose_from_file(str(full_name))
df = pd.read_csv(meta_file)
n_resi, c_resi = df["n_resi"].values[0], df["c_resi"].values[0]
if c_resi > 0:
    pose.delete_residue_range_slow(pose.size() - c_resi + 1, pose.size())
if n_resi > 0:  # else keep the N termini
    pose.delete_residue_range_slow(1, n_resi)
full_name = native_relaxed_path.parent.parent.joinpath(f"native/{{pdbid}}_cut.pdb")
pose.dump_pdb(str(full_name))""")

python = "/home/mads/miniconda3/envs/shapedesign/bin/python"
def setup_from_af(experiment_number, af_data_path, production_path, structures, slide, maxiter, popsize,
                  mutate, recombination, flex, low_memory, normalize_score, scheme, seed,
                  local_search,  number_of_runs, wall_time, redo_only_configs, selection, direction, allow_flip, init_input_fix_percent):

    exp_path = Path(production_path).joinpath(experiment_number)
    exp_path.mkdir(exist_ok=True, parents=True)
    paths = [path for path in Path(af_data_path).glob("**/*.csv") if path.stem in structures and all(part not in path.name for part in ("meta", "xtrans", "pack_info"))]
    # all_scripts = []
    dfdir = pd.read_csv(direction)
    for path in paths:
        df_meta = pd.read_csv(path.parent.joinpath(path.stem + "_meta.csv"))
        df_all = pd.read_csv(path)
        pdbid = path.stem
        # --- create paths
        pdb_path = exp_path.joinpath(pdbid)
        pdb_path.mkdir(exist_ok=True)
        native_path = pdb_path.joinpath("native")
        native_path.mkdir(exist_ok=True)
        input_path = pdb_path.joinpath("inputs")
        input_path.mkdir(exist_ok=True)
        subunits_path = pdb_path.joinpath("subunits")
        subunits_path.mkdir(exist_ok=True)
        # preperation_path = pdb_path.joinpath("preperation")
        # preperation_path.mkdir(exist_ok=True)
        run_dir = pdb_path.joinpath("run")
        run_dir.mkdir(exist_ok=True)
        config_dir = run_dir.joinpath("configs")
        config_dir.mkdir(exist_ok=True)
        results_dir = run_dir.joinpath("results")
        results_dir.mkdir(exist_ok=True)
        log_dir = run_dir.joinpath("logs")
        log_dir.mkdir(exist_ok=True)

        # transfer the subunits to the subunits folder
        direction = dfdir[dfdir["pdb"] == pdbid]["direction"].values[0]
        for file in path.parent.parent.glob(f"pdbs/{direction}/*"):
            file_to = subunits_path.joinpath(file.name)
            shutil.copy(str(file), str(file_to))

        # get the native file locations
        native_file_src = list(Path(SYMMETRICAL).glob(f"**/crystal_repr/native/{pdbid}_crystal.pdb"))[0]
        native_file_symmetric_src = list(Path(SYMMETRICAL).glob(f"**/input/native/{pdbid}.cif"))[0]
        native_symdef_src = list(Path(SYMMETRICAL).glob(f"**/symdef/native/{pdbid}.symm"))[0]

        # Read in the posese and cut the structures as they are in the predictions
        init("-initialize_rigid_body_dofs true")
        native_file = pose_from_file(str(native_file_src))
        native_file_symmetric = pose_from_file(str(native_file_symmetric_src))
        n_resi = df_meta["n_resi"].values[0]
        c_resi = df_meta["c_resi"].values[0]
        native_file = AF2EvoDOCK.cut_multimer_poses(native_file, n_resi=n_resi, c_resi=c_resi)
        native_file_symmetric = AF2EvoDOCK.cut_monomeric_pose(native_file_symmetric, n_resi=n_resi, c_resi=c_resi)

        # put them into their new locations
        native_file_new_loc = native_path.joinpath(native_file_src.name)
        native_file_symmetric_new_loc = native_path.joinpath(native_file_symmetric_src.name)
        native_symdef_new_loc = native_path.joinpath(native_symdef_src.name)

        native_file.dump_pdb(str(native_file_new_loc))
        native_file_symmetric.dump_pdb(str(native_file_symmetric_new_loc))
        # just copy the symdef directly
        shutil.copy(native_symdef_src, native_symdef_new_loc)

        # we need to put the normalized symmetry file into the input folder
        symmetry = native_file_src._str.replace("/home/shared/databases/SYMMETRICAL/", "")[0] # insane hack to get the symmetry file
        cn = str(df_all["mer"].values[0]) #
        input_symdef_src = AF2EvoDOCK.cn_to_sympath(symmetry, cn)
        input_symdef_new_loc = input_path.joinpath(input_symdef_src.name)
        shutil.copy(input_symdef_src, input_symdef_new_loc)



        # NATIVE RELAX
        # with open(preperation_path.joinpath("relax_native_script.sh"), "w") as f:
        #     f.write("#! /bin/bash\n")
        #     symmetric_relax = "/home/mads/miniconda3/envs/shapedesign/bin/symmetric_relax.py"
        #     # shutil.copy(native_pdb, native_relaxed_name)
        #     f.write(f"cd {native_relaxed_path} && {python} {symmetric_relax} --file {native_file_symmetric} --symmetry_file {native_symdef}")
        #     f.write(f" --cycles 15 --ex1 --ex2aro --constrain_coords --bb --sc --jumps --move_jumpdofs JUMPHFfold1:z:2000 JUMPHFfold111:x:2000\n")
        #     f.write(f"{python} {preperation_path.joinpath('cut_script.py')}")
        #     # # run the
        #     #

        # SUBUNITS RELAX, BACKRUB, NMA
        # figure out how many to run in total
        # i = 0
        # runs = {"relax_w_constrains": [], "relax": [], "nma": [], "backrub": []}
        # # first add all AlphaFold instsances to relax_w_constraints
        # for model_name_cut in df.sort_values("avg_plddt")["model_name_cut"]:
        #     if not i < ensemble_size:
        #         break
        #     else:
        #         runs["relax_w_constrains"].append(model_name_cut)
        #         i += 1
        # while i < ensemble_size:
        #     # Then add unitl
        #     for model_name_cut in df.sort_values("avg_plddt")["model_name_cut"]:
        #         if not i < ensemble_size:
        #             break
        #         else:
        #             sample_type = np.random.choice(["relax", "nma", "backrub"], p=[0.30, 0.40, 0.30])
        #             runs[sample_type].append(model_name_cut)
        #             i += 1
        # # create all the python scripts
        # if not redo_only_configs:
        #     all_scripts += create_flexbb_creater_script(preperation_path, subunits_path, native_file_symmetric_cut, runs, ensemble_size)
        # # here you have to implement what is in  run

        # create copy script
        # make_cut_and_copy_script(native_file_symmetric, pdbid, meta_file=meta)

        # bounds are set based on the symmetry and the CN
        if symmetry == "T":
            if cn == "2":
                angle_z = 90
            elif cn == "3" or cn == "HF":
                angle_z = 60
            else:
                raise ValueError
        elif symmetry == "O":
            if cn == "2":
                angle_z = 90
            elif cn == "3":
                angle_z = 60
            elif cn == "4" or cn == "HF":
                angle_z = 45
            else:
                raise ValueError
        elif symmetry == "I":
            if cn == "2":
                angle_z = 90
            elif cn == "3":
                angle_z = 60
            elif cn == "5" or cn == "HF":
                angle_z = 36
            else:
                raise ValueError
        bounds = f"1000,{angle_z},5,40,40,40"
        init_bounds = f"0,{angle_z},5,40,40,40"

        # this has tobe run one more time in order to create this
        rel = Path("../")
        ini_files = []
        for i in range(number_of_runs):
            # create sesult dir:
            ini_file = f"{pdbid}{i}.ini"
            ini_files.append(ini_file)
            results_dir_i = results_dir.joinpath(str(i))
            results_dir_i.mkdir(exist_ok=True)
            create_from_flexbb_template(native_input=rel.joinpath(*native_file_new_loc.parts[-2:]),
                                        native_symmetric_file=rel.joinpath(*native_file_symmetric_new_loc.parts[-2:]),
                                        native_symdef_file=rel.joinpath(*native_symdef_new_loc.parts[-2:]),
                                        input_symdef_path=rel.joinpath(*input_symdef_new_loc.parts[-2:]),
                                        output_path=rel.joinpath(f"results/{i}/"),
                                        config_name=config_dir.joinpath(ini_file),
                                        subunits=rel.joinpath("subunits"),
                                        scheme=scheme,
                                        bounds=bounds,
                                        slide=slide,
                                        low_memory=low_memory,
                                        normalize_score=normalize_score,
                                        popsize=popsize,
                                        mutate=mutate,
                                        recombination=recombination,
                                        maxiter=maxiter,
                                        selection=selection,
                                        init=init_bounds,
                                        allow_flip=allow_flip,
                                        init_input_fix_percent=init_input_fix_percent)
        jobname = run_dir.joinpath(f"{pdbid}_job.sh")
        create_evodock_submission_script_lunarc(jobname, ini_files, jobname=f"{pdbid}", time=wall_time)

    # if not redo_only_configs:
    #     run_path = exp_path.joinpath("run")
    #     run_path.mkdir(exist_ok=True)
    #     nodefile_path = run_path.joinpath("nodefile.txt")
    #     create_node_file(nodefile_path, mailer=34)
    #     parallel_input_file = run_path.joinpath("input_file.txt")
    #     create_parallel_input_file(parallel_input_file, all_scripts)
    #     parallel_script_path = run_path.joinpath("ensemble_prep.sh")
    #     command = f"{python} {{1}}"
    #     joblog_path = run_path.joinpath("joblog")
    #     create_parallel_script(parallel_script_path, command, nodefile_path, joblog_path, parallel_input_file)

def main():

    description = "Sets up EvoDOCK runs from the AF ensemble"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--experiment_number', help="Number of the experiment", type=str, required=True)
    parser.add_argument("--production_path", help="path to where to produce all the outputs.", type=str, required=True)
    parser.add_argument("--af_data_path", help="path to the af csv files that was generated to an ensemble", type=str, required=True)
    parser.add_argument("--number_of_runs", help="number of runs to do", type=int, required=True)
    parser.add_argument("--structures", help="Which structures to use", nargs="+", type=str, required=True)
    # parser.add_argument("--ensemble_size", help="size of the ensemble", type=int, default=150)
    parser.add_argument("--slide", help="Use slide in EvoDOCK", default=True, type=lambda x: bool(strtobool(x)), nargs='?')
    parser.add_argument("--maxiter", help="", type=int, default=100, required=False)
    parser.add_argument("--popsize", help="", type=int, default=100, required=False)
    parser.add_argument("--mutate", help="", type=float, default=0.1, required=False)
    parser.add_argument("--recombination", help="", type=float, default=0.7, required=False)
    parser.add_argument("--flex", help="To use flexible docking or not", action="store_true")
    parser.add_argument("--low_memory", help="To use low memory or not", action="store_true")
    parser.add_argument("--normalize_score", help="To normalize the score or not", default=False, type=lambda x: bool(strtobool(x)), nargs='?')
    parser.add_argument("--selection", help="selection type to use", type=str, default="total", choices=["interface", "total"])
    parser.add_argument("--scheme", help="scheme type to use", type=str, default="RANDOM", choices=["BEST", "RANDOM"])
    parser.add_argument("--direction", help="direction to use", type=str, required=True)
    #
    parser.add_argument("--allow_flip", help="allow flipping or not", type=lambda x: bool(strtobool(x)), nargs='?', default=False)
    parser.add_argument("--init_input_fix_percent", help="The amount of structures in percent where the 4 dofs you get from alphafold "
                                                         "will not be perturbed.", type=float, default=0)
    # init_input_fix_percent
    # allow_flip
    # FIXME: make random instead!
    parser.add_argument("--seed", help="Sets a constant seed and ofsets the individual runs by 12345", type=int)
    # parser.add_argument("--dope_library", help="Use the doped Input in the BB library for FlexBB",
    #                     default=False, type=lambda x: bool(strtobool(x)), nargs='?')
    # parser.add_argument("--rmsd", help="The flexbb subunits cannot have less than the set RMSD to the native structure", type=float)
    parser.add_argument("--local_search", help="local search type to use", type=str, default="symshapedock")
    # parser.add_argument("--rmsd_acceptance", help="local search type to use", type=float, default=1.0)
    # parser.add_argument("--dssp_acceptance", help="local search type to use", type=float, default=5.0)
    parser.add_argument("--wall_time", help="wall time", type=str)
    parser.add_argument("--redo_only_configs", help="only changes the configs", action="store_true")
    args = parser.parse_args()

    setup_from_af(args.experiment_number, args.af_data_path, args.production_path, args.structures, args.slide, args.maxiter, args.popsize,
                  args.mutate, args.recombination, args.flex, args.low_memory, args.normalize_score, args.scheme, args.seed,
                  args.local_search,  args.number_of_runs, args.wall_time, args.redo_only_configs, args.selection, args.direction,
                  args.allow_flip, args.init_input_fix_percent)

if __name__ == "__main__":
    main()