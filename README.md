
# EvoDOCK 

EvoDOCK is a software for Heterodimeric and Symmetric Protein-Protein docking.

Heterodimeric docking has been published at:
[A memetic algorithm enables global all-atom protein-protein docking with sidechain flexibility](https://www.cell.com/structure/fulltext/S0969-2126(22)00372-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0969212622003720%3Fshowall%3Dtrue)

Symmetric docking has been published (as preprint) at:
[Accurate prediction of protein assembly structure by combining AlphaFold and symmetrical docking](https://www.biorxiv.org/content/10.1101/2023.06.22.546069v1)

# System Requirements

## Hardware requirements

EvoDOCK only requires a standard computer to run but as population size increase more RAM is needed. 

## Software requirements

### OS Requirements

This package is supported for Linux (macOS will soon be added). The package has been tested on the following systems:

Linux: Ubuntu 20.04.5-6 and CentOS Linux 7

## Package requirements

If using heterodimeric only the following packages must be installed: 
* Python-3.6 or later (PyRosetta dependency). 
* PyRosetta (http://www.pyrosetta.org) (Can be installed with Anaconda). You need to obtain a license before use (see the link) 

If using Symmetric Protein-Protein docking these packages must be installed:
* MAFFT (https://mafft.cbrc.jp/alignment/software/) (can be installed with Anaconda/brew/apt)
* mpi4py and its requirements (https://mpi4py.readthedocs.io/en/stable/install.html) (can be install with Anaconda/pip)
  
Furthermore the following packages are also needed but are automatically installed by pip install using the setup.py script (see Installation Guide).

```
cubicsym @ git+https://github.com/Andre-lab/cubicsym,
cloudcontactscore @ git+https://github.com/Andre-lab/cloudcontactscore,
numpy>=1.21.0
pandas>=1.3.4
pillow>=9.1.0
scipy>=1.7.1
seaborn>=0.11.2
setuptools>=44.0.0
imageio>=2.10.1
matplotlib>=3.4.3
```

# Installation Guide

Clone the cubicsym repository and `cd` into it. Then run the install script.
```console
git clone https://github.com/Andre-lab/evodock.git
cd ./evodock
pip install .
```

Installation time should only take a couple of seconds but downloading the required pacakges and installing them can take several minutes.

# Usage

## Preprocessing a complex with prepacking
Before running EvoDOCK it is important to prepack the input files as so (takes several seconds): 

```console
python ./scripts/prepacking.py --file <input_file>
```

## Converting AlphaFold predictions to a symmetric EvoDOCK ensemble
 
`scripts/af_to_evodock.py` converts AlphaFold 2 and AlphaFold-Multimer predictions to an EvoDOCK ensemble.
The script is well documented. Use `python scripts/af_to_evodock.py -h` to see more. The output will already be prepacked.

Below, 2 examples of running the script for creating an ensemble for Reassembly docking or Complete assembly docking is given. You need to download `af_data.tar.gz` [here](https://zenodo.org/record/8047514). Unzip it with 

```console
tar -xf af_data.tar.gz
```

Put the AF_data in `evodock/inputs` before running the tests below. 

Preparing an ensemble for Reassembly docking (takes a few minutes):
```console
scripts/af_to_evodock.py --path inputs/AF_data/local --symmetry O --ensemble Local --out_dir tests/outputs/ --max_multimers 5 --max_monomers 5 --modify_rmsd_to_reach_min_models 50 --max_total_models 5 --align_structure inputs/input_pdb/3N1I/3N1I.cif 
```

Preparing an ensemble for Complete assembly docking (takes a few minutes):
```console
scripts/af_to_evodock.py --path inputs/AF_data/globalfrommultimer --symmetry T --ensemble GlobalFromMultimer --out_dir tests/outputs/ --max_multimers 5 --max_monomers 5 --modify_rmsd_to_reach_min_models 50 --max_total_models 5
```

## EvoDOCK 
EvoDOCK can be run with different configurations given a specifc config.ini input file as so: 

```console
python evodock.py configs.ini
```

The following sections describe how to configure EvoDOCK through the config file. These options are available: 
1. [Docking]
2. [Input]
3. [Outputs]
4. [DE]
5. [Flexbb] 
6. [Bounds]
7. [Pymol]
8. [RosettaOptions]
9. [Native]

Examples of config files for different EvoDOCK configurations are found in the `config` folder with the following behavior: 

1. Heteromeric docking with single ligand and receptor backone (takes a few minutes): `configs/heterodimeric/sample_dock_single.ini` 
2. Heteromeric docking with flexible backbones (takes a few minutes): `configs/heterodimeric/sample_dock_flexbb.ini` 
3. Reassembly docking with a single backbone (takes a few minutes): `configs/symmetric/reassembly_single.ini` 
4. Reassembly docking with flexible backbones (takes a few minutes): `configs/symmetric/reassembly_flexbb.ini` 
5. Complete assembly docking with flexible backbones (takes a few minutes): `configs/symmetric/complete_assembly.ini`  

### 1. [Docking]
Specifies the type of docking protocol used of which there are 3 options:
1. `Local` For heterodimeric local docking AND symmetric Reassembly docking
3. `Global` For heterodimeric global docking 
4. `GlobalFromMultimer` For symmetric Complete assembly docking
```dosini
[Docking]
type=<Local/Global/GlobalFromMultimer>
```

### 2. [Input]

Specifies the input type. For heteromeric docking you need to specify the either
single or `ligands` AND `receptors` for docking either 2 single backbones or 2 sets of multiple backbones. For heterodimeric
docking a `template` can be supplied. This is used to extact rotamers and to initially align the receptor and ligand onto.
```dosini
[Input]
single=<path to a pdb file containing containg the heterodimer (2 chains)>
```
or
```dosini
[Input]
ligands=<path to ligands>
receptors=<path to receptors>
```
or
```dosini
[Input]
template=<path to a pdb file to serve as a template>
ligands=<path to ligands>
receptors=<path to receptors>
```

For symmetric docking you need to specify the `symdef_file` and either `single` or `subunits` for docking either a single or multiple backbones
```dosini
[Input]
single=<path to single pdb file>
symdef_file=<path to the symdef file>
```
or 
```dosini
[Input]
subunits=<path to a directory containing all of the subunits>
symdef_file=<path to the symdef file>
```

### 3. [Outputs]

Output options for the results:
1. `output_path` Directory in which to output all files.
2. `output_pdb` Output pdb
3. `output_pdb_per_generation` Output the best pdb for each generation

```dosini
[Outputs]
output_path=<path to the output directory>
output_pdb=<boolean>
output_pdb_per_generation=<boolean>
```

#### 4. [DE]
Differential Evolution options:
1. `scheme`: The selection strategy for the base vector at mutation operation. 1. Selection randomly (=RANDOM, default), 2. Select the best (=BEST)
2. `popsize`: The size of the population. Default is 100.
3. `mutate`: mutation rate (weight factor F). Must be between 0 and 1.0. Default is 0.1.
4. `recombination`: crossover probability (CR). Must be between 0 and 1.0. Default is 0.7.
5. `maxiter`: Generations to perform. Default is 50
6. `local_search`: The local search docking scheme. For heteromeric docking use [None, only_slide, mcm_rosetta] for symmetryic docking use symshapedock. Default for heterodimeric docking is mcm_rosetta and for symmetryic docking symshapedock.
7. `slide`: Use sliding or not. Default is True.
8. `selection`: The energy type for use when selecting. 1. Select by interface (=interface, default for symmetric docking), 2. select by total energy (=total, default for heterodimeric docking)

```dosini
[DE]
scheme=<RANDOM/BEST>
popsize=<integer>
mutate=<float>
recombination=<float>
maxiter=<integer>
local_search=<None/only_slide/mcm_rosetta/symshapedock>
slide=<boolean>
selection=<interface/total>
```

### 5. [Flexbb]

If this section is present EvoDOCK will do flexible backbone docking. 2 options can be set:
1. `swap_prob` The probability of doing a backbone trial. Must be in the interval: [0, 1.0] 
2. `low_memory_mode` Will save memory by only loading in 1 backbone at the time at the cost of some computional time. Is only available for symmetrical docking and is highly recommend when using symmetrical docking. 

```dosini
[Flexbb]
swap_prob=<float>
low_memory_mode=<boolean>
```

### 6. [Bounds]
Set options for the bounds of the rigid body parameters when doing symmetrical docking:
1. `init`: The initial bounds the rigid body parameters are sampled in; [z, λ, x, ψ, ϴ, φ] for cubic symmetric docking.
2. `bounds:`:  The maximum bounds the rigid body parameters are sampled in; [z, λ, x, ψ, ϴ, φ] for cubic symmetric docking.
3. `init_input_fix_percent`: The percent chance of keeping an individual to its exact input values and not randomizing inside the init bounds. Should be between 0 and 100. 
4. `allow_flip`: allow the individual to flip 180 degrees. 
5. `xtrans_file`: The path to the file containing the x translations for each subunit. This file is output from af_to_evodock.py when running with --ensemble=GlobalFromMultimer

```dosini
[Bounds]
init=<initial bounds, example: 0,60,5,40,40,40>
bounds=<maximum bounds, example: 1000,60,5,180,40,40>
init_input_fix_percent=<float>
allow_flip=<boolean>
xtrans_file=<path to the xtrans_file>
```

### 7. [Pymol]
EvoDOCK can be run with PyMOL as described in https://www.rosettacommons.org/docs/latest/rosetta_basics/PyMOL.
This sets options for PyMOL:
1. `on`: Use PyMOL.
2. `history`: Turn history on.
3. `show_local_search`: Show the processes in during local search.
4. `ipaddress`: The IP address to use.  

```
[Pymol]
on=<boolean>
history=<boolean>
show_local_search=<boolean>
ipaddress=<IP address>
```

### 8. [RosettaOptions]
Rosetta flags to use. Any can be specified. 

```dosini
[RosettaOptions]
initialize_rigid_body_dofs=<boolean>
```

### 9. [Native]

Calculates metrics againts the native structure (RMSD for instance). There are 3 input types:
1. `crystallic_native` The native structure
2. `symmetric_input` The symmetric input file of the native structure
3. `symdef_file` The symdef file for the native structure
4. `lower_diversity_limit` The lowest RMSD limit the structures should have to their native structure

2 and 3 is required for symmetric docking.

```dosini
[Native]
crystallic_input=<path to native structure>
symmetric_input=<path to the symmetric input>
symdef_file=<path to the input file>
lower_diversity_limit=<float>
```

# Interpretation of the EvoDOCK output:

EvoDOCK produces several different log files:

1. `evolution.csv` is a general summary of the evolutionary process across the entire population. Per generation (gen) it logs:
    - The average energy of the population (avg)
    - The lowest energy of population (best)
    - The RMSD of the best individual with the lowest energy (rmsd_from_best) if running with a native structure.

2. `popul.csv` is a general summary of the evolutionary process for each individual in the population. Per generation (gen) it logs:
    - The current score (sc_*)
    - The current rmsd (rmsd_*) if running with a native structure.
    - The current interface score (Isc_*)
    - The current Interface rmsd (Irmsd_*) if running with a native structure.

3. `trials.csv` is the equivalent file to popul.csv, but it reports the trials (candidates) generated during the each generation. This can be practically useful in case that you want to check if the DE+MC is creating proper candidates that can contribute to the evolution.

4. `time.csv` is the computational time (in seconds) for each generation (gen).

5. `all_individual.csv` contains, for each generation (gen), the best genotype (rigid body degrees of freedom) of all individuals. 

6. `best_individual.csv` contains, for each generation (gen), the best genotype (rigid body degrees of freedom) of the individual with lowest energy value.

7. `population_swap.csv` contains, for each generation (gen), the backbone swap success. 

8. `flip_fix.csv` list for each individual if they were initially flipped or fixed. Is useful when running running with `GlobalFromMultimer` .

9. `ensemble.csv` contains, for each generation (gen), the name of file that is used as the current backbone for each individual.

EvoDOCK also outputs a pdb file of the final optimized model (`final_docked_evo.pdb`) as well as the lowest energy structure for each geneation (`evolved.pdb`)

# Symmetric relax of EvoDOCK output structures

The script `scripts/symmetric_relax.py` can be used to relax structures from the EvoDOCK output. The script is well documented. Use `python scripts/symmetric_relax.py -h` to see more.
It is advisable to use this script when parsing AlphaFold models as compared to Rosettas relax protocol as it guards against the structures blowing up if the AlphaFold structures have bad energies. It does however require more user interference as explained below.

Use the `ensemble.csv` and `all_individual.csv` to get the backbone and genotype ([z, λ, x, ψ, ϴ, φ]) for any model you want for any generation. If you want the lowest energy one parse this into a pandas 
DataFrame and fish it out. Use the genotype to modify the symmetry file you parsed to EvoDOCK (see https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry for more information about the symmetry files in Rosetta). 
You have to modify the set_dof lines to match the genotype. 
Then to the script parse the backbone (pdb file) to --file and the modified symmetry file to --symmetry_file.

A test can be run with (can take up to an hour or more):

```console
python scripts/symmetric_relax.py --file inputs/input_pdb/2CC9/2CC9_tobe_relaxed.pdb --symmetry_file inputs/symmetry_files/2CC9_tobe_relaxed.symm --rosetta_out tests/outputs/ --input_out tests/outputs/ 
```

# Differential Evolution Algorithm

Differential Evolution [Price97] is a population-based search method. DE creates new candidate solutions by combining existing ones according to a simple formula of vector crossover and mutation, and then keeping whichever candidate solution has the best score or fitness on the optimization problem at hand.

# Bibliography

* Storn, R., Price, K. Differential Evolution – A Simple and Efficient Heuristic for global Optimization over Continuous Spaces. Journal of Global Optimization 11, 341–359 (1997). https://doi.org/10.1023/A:1008202821328 
