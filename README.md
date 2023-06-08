
# EvoDOCK 

EvoDOCK is a software for Heterodimeric and Symmetric Protein-Protein docking.

Heterodimeric docking has been published at:
[A memetic algorithm enables global all-atom protein-protein docking with sidechain flexibility](https://www.biorxiv.org/content/10.1101/2021.04.12.437963v3)

Symmetric docking has been published at:
[missing](link:missing)

# Installation 

clone the evodock repository and ```cd``` into it
```
git clone https://github.com/Andre-lab/evodock.git
cd ./evodock
```

## Conda installation
EvoDOCK can be installed with [Anaconda](https://docs.anaconda.com/free/anaconda/install/index.html).

A license for PyRosetta is needed and can be obtained from https://els2.comotion.uw.edu/product/pyrosetta. 
Use the provided USERNAME and PASSWORD and insert it into the following line in ```env.yml```:
```
- https://USERNAME:PASSWORD@conda.rosettacommons.org
```

Then run the following conda command (it will take a few minutes):

```
conda install evodock
```

## pip installation

This package is only compatible with Python-3.6 or later (PyRosetta dependency)

* Download and install PyRosetta http://www.pyrosetta.org/dow
* Download and install MAFFT (If using Symmetry) (https://mafft.cbrc.jp/alignment/software/)
* Install the package itself:

```console
git clone https://github.com/Andre-lab/evodock.git
cd evodock
pip install -r requirements.txt
```

or

```console
git clone https://github.com/Andre-lab/evodock.git
pip setup.py install 
```

or 

```console
pip install git+https://github.com/Andre-lab/evodock.git
```

A setup.py and environment.yml files are provided to use alternative installation using pip or conda.


# Basic Usage

0. Preprocess complex pdb with prepacking

```console
python ./scripts/prepacking.py <input_pdb>
```


1. Create a configuration file following the example found at [sample\_dock.ini](https://github.com/Andre-lab/evodock/blob/2fbc755cf84f64641153ad75757ad4bb3bf6ff3f/configs/sample_dock.ini)

```dosini
[Docking]
# selects docking protocl [Global, Local]
type=Global

[Inputs]
# complex pdb
pose_input=/inputs/input_pdb/1ACB/1ACB_c_u_0001.pdb
native_input=/inputs/native_pdb/1ACB/1ACB_c_b.pdb

[Outputs]
# output file log
output_path=sample_dock/
output_pdb=True

[DE]
# evolution algorithm parent strategy [RANDOM, BEST] 
scheme=BEST
# population size
popsize=10
# mutation rate (weight factor F) 
mutate=0.9
# crossover probability (CR) 
recombination=0.3
# maximum number of generations/iterations (stopping criteria)
maxiter=10
# hybrid local search strategy [None, only_slide, mcm_rosetta]
local_search=mcm_rosetta

```
information about the DE parameters can be found at https://en.wikipedia.org/wiki/Differential_evolution


2. Run with the algorithm with the desired configuration

```console
python evodock.py configs/sample_dock_global.ini
```

or 

```console
python -m evodock configs/sample_dock_global.ini
```


## Configuration Details

Files configs/sample\_dock\_global.ini, configs/sample\_dock\_flexbb.ini and configs/sample\_dock_refinement.ini contains configuration examples for Global Docking, Flexible Backbone Docking and Global Docking with and initial population. There's also an configs/sample\_sym_dock.ini that contains a configuration example for running EvoDock with symmetry. For more information see **Running with symmetry** further below.  

### Section [Inputs]

At pose\_input, you might provide the path to a complex with two chains, which previously was preprocessed with a prepack protocol in order to fix possible collisions at the sidechain. An script at script folders is provided. 

### Section [Outputs]

At  output\_path indicates the output folder for the results in .csv format.
output\_pdb is a boolean to dump pdbs during the evolution and the final evolved protein.


### Section [Docking]
Option "type" allows to select between global docking (Global), local docking (Local), flexible backbone (Flexbb) and
using an starting population such as models from ClusPro (Refinement).

### Section [DE]
The set of parameters for Differential Evolution ([DE])  that you must change for a production run are populsize (from 10 to 100) and maxiter (from 10 to 100), which would lead into an evolution of 100 individuals during 100 iterations/generations. Evolutionary parameters (mutation F and crossover CR), can be fine tuned for specific purposes, although this set (0.3 and 0.9) have shown a good balance between exploration and exploration at our benchmark runs, which leads into good results. Scheme corresponds to the selection strategy for the base vector at mutation operation (https://en.wikipedia.org/wiki/Differential_evolution for more details). Parameter "local\_search" can be changed to None (aka, only DE is performed), only\_slide (local search operation is equivalent to apply slide\_into\_contact) or mcm\_rosetta (which applies slide\_into\_contact + MC energy minimization and sidechain optimization, recommended option and used at our benchmarks)

### Section [Flexbb] (optional for Docking type "Flexbb")
Uses path\_ligands and path\_receptors to indicate the path of *.pdb files with different backbone ensembles.

### Section [Refine] (optional for Docking type "Refine")
Uses init\_pdbs to indicate the path of *.pdbs used as initial population, i.e. models from ClusPro.

## Running with symmetry
To run with symmetry add the following configuration as below.

```dosini
[Symmetry]
input_symdef_file=/inputs/symmetry_files/1stm.symm
native_symdef_file=/inputs/symmetry_files/1stm.symm
symdofs=JUMP5fold1:z:angle_z,JUMP5fold111:x,JUMP5fold1111:angle_x:angle_y:angle_z
symbounds=15,15,15,15,15,15
initialize_rigid_body_dofs=true
```


# Interpret output:

It is going to produce 4 different log files:

-   evolution\*csv is a summary of the evolutionary process, which indicates the number of generation,

average energy of the population, lowest energy of population and the RMSD of the best individual with the lowest energy.

-   popul\*csv is the status of each generation during the evolution. Each line correponds to the population information of one generation.

-   interface\*csv is similar to popul, but it reports the interface energy value and the iRMSD for each corresponding individual at each generation.

-   trials\*csv is the equivalent file to popul\*csv, but it reports the trials (candidates) generated during the each generation. This can be practically useful in case that you want to check if the DE+MC is creating proper candidates that can contribute to the evolution.

-   time\*csv is the computational time (in seconds) for each generation.

-   best\*csv contains, at each line, the rotation (first 3 values) and translation (3 values) of the individual with lowest energy value.


## Getting images

### Get scatter plot

python ./scripts/make\_scatter\_plot.py "<path\_to\_popul\*.csv>"

It creates the global energy value vs RMSD plot if input is popul*csv or interface energy vs iRMSD plot if input corresponds to interface*csv. Each point corresponds to an individual in the last generation. Several \*csv files can be specified in order to collect the results from different independent runs, where each color corresponds to a run.

![interface Energy vs iRMSD scatter plot](https://github.com/Andre-lab/evodock/blob/main/images/scatterplot.png)


### Get evolution performance

For each popul\*csv

python ./scripts/make\_evolution\_plot.py <path to evolution\*.csv>

Creates a lineplot where y-axis corresponds to the global energy function (used as fitness function during the evolution) and x-axis corresponds to each generation.

Green line corresponds to the average energy value of the population, while the red line corresponds to the lowest energy value of the population. A proper evolution should maintain a close distance between both lines and average line should follow the tend of the lowest energy line. That would indicate that the population evolves towards the best energy individual. In case that there is a large different between both lines, F and CR parameters should be tuned. For example, reducing the exploration of the algorithm by decreasing the value of F.

![Evolution sample image](https://github.com/Andre-lab/evodock/blob/main/images/evolution_sample.png)

# Differential Evolution Algorithm

Differential Evolution [Price97] is a population-based search method. DE creates new candidate solutions by combining existing ones according to a simple formula of vector crossover and mutation, and then keeping whichever candidate solution has the best score or fitness on the optimization problem at hand.


# Bibliography

* Storn, R., Price, K. Differential Evolution – A Simple and Efficient Heuristic for global Optimization over Continuous Spaces. Journal of Global Optimization 11, 341–359 (1997). https://doi.org/10.1023/A:1008202821328 
