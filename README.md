
# Protein-Protein docking using a Memetic Algorithm: EvoDOCK

:warning: under construction


Repository corresponding to the code used at article: 

[A memetic algorithm enables global all-atom protein-protein docking with sidechain flexibility](https://www.biorxiv.org/content/10.1101/2021.04.12.437963v1)

# Dependencies

* PyRosetta==4
* numpy==1.19.4



# Installation

This package is only compatible with Python 3.4 and above. To install this package, please follow the instructions below:

* Install the previous descripted dependencies
* Download and install PyRosetta following the instructions found at http://www.pyrosetta.org/dow
* Install the package itself:

```console
git clone https://github.com/Andre-lab/evodock.git
cd evodock
pip install -r requirements.txt
```

# Basic Usage

0. Preprocess complex pdb with prepacking

```console
python ./scripts/preprocessing <input_pdb>
```


1. Create a configuration file following the example found at [sample\_dock.ini](https://github.com/Andre-lab/evodock/blob/2fbc755cf84f64641153ad75757ad4bb3bf6ff3f/configs/sample_dock.ini)

```dosini
[inputs]
# complex pdb
pose_input=./inputs/pdbs/1ppe_IE.prepack.pdb

[outputs]
# output file log
output_file=quick_evolution_sample.log

[DE]
# evolution algorithm parent strategy [RANDOM, BEST] 
scheme=BEST
# population size
popsize=10
# mutation rate (weight factor F) 
mutate=0.9
# crossover probability (CR) 
recombination=0.6
# maximum number of generations/iterations (stopping criteria)
maxiter=10
# hybrid local search strategy [None, only_slide, mcm_rosetta]
local_search=mcm_rosetta

```
information about the DE parameters can be found at https://en.wikipedia.org/wiki/Differential_evolution


2. Run with the algorithm with the desired configuration

```console
python evodock.py configs/sample_dock.ini
```


## Configuration Details

### Section [inputs]

At pose<sub>input</sub>, you might provide the path to a complex with two
chains, which previously was preprocessed with a prepack protocol in order to fix possible collisions at the sidechain. An
script at script folders is provided. 

### Section [outputs]

At  output<sub>file</sub> you can indicate the output folder and log file for
the output. Log file name should contain the word "evolution" to not produce
errors (todo: improve the code of this). The name of the output<sub>file</sub> should be
different for each independent run. 

### Section [DE]
The set of parameters for Differential Evolution ([DE])  that you must change for a production run are populsize (from 10 to 100) and maxiter (from 10 to 100),
which would lead into an evolution of 100 individuals during 100
iterations/generations. Evolutionary parameters (mutation F and crossover CR),
can be fine tuned for specific purposes, although this set (0.3 and 0.9) have
shown a good balance between exploration and exploration at our benchmark runs,
which leads into good results. Scheme corresponds to the selection strategy for
the base vector at mutation operation (https://en.wikipedia.org/wiki/Differential_evolution for more details).
Parameter "local<sub>search</sub>" can be changed to None (aka, only DE is performed),
only<sub>slide</sub> (local search operation is equivalent to apply slide<sub>into</sub><sub>contact</sub>)
or mcm<sub>rosetta</sub> (which applies slide<sub>into</sub><sub>contact</sub> + MC energy minimization and
sidechain optimization, recommended option and used at our benchmarks)


# Interpret output:

It is going to produce 4 different log files:

-   evolution\*log is a summary of the evolutionary process, which indicates the number of generation,

average energy of the population, lowest energy of population and the RMSD of the
best individual with the lowest energy.

-   popul\*log is the status of each generation during the evolution. It contains

two lines per generation: one corresponding to the energy value of each
individual and other to the corresponding RMSD.

-   interface\*log is similar to popul, but it reports the interface energy value

and the iRMSD for each corresponding individual at each generation.

-   trials\*log is the equivalent file to popul\*log, but it reports the trials

(candidates) generated during the each generation. This can be practically
useful in case that you want to check if the DE+MC is creating proper
candidates that can contribute to the evolution.

-   time\*log is the computational time (in seconds) for each generation.

-   best\*log contains, at each line, the rotation (first 3 values) and

translation (last 3 values) of the individual with lowest energy value. 


## Getting images

### Get scatter plot

python ./scripts/make<sub>scatter</sub><sub>plot.py</sub>
<path<sub>to</sub><sub>popul</sub>\*.log>

It creates the interface energy vs iRMSD plot,
commonly used to evaluate the sampling results. Each point
corresponds to an individual in the last individual. Several popul\*log files
can be specified in order to collect the results from different independent runs.


### Get evolution performance

For each popul\*log

python ./scripts/make<sub>evolution</sub><sub>plot.py</sub> <path to evolution\*.log>

Creates a lineplot where y-axis corresponds to the global energy function (used
as fitness function during the evolution) and x-axis corresponds to each
generation.

Green line corresponds to the average energy value of the population, while the
red line corresponds to the lowest energy value of the population. A proper
evolution should maintain a close distance between both lines and average line
should follow the tend of the lowest energy line. That would indicate that the
population evolves towards the best energy individual. In case that there is a
large different between both lines, F and CR parameters should be tuned. For
example, reducing the exploration of the algorithm by decreasing the value of F.



# Differential Evolution Algorithm

Differential Evolution [Price97] is a population-based search method. DE creates new candidate solutions by combining existing ones according to a simple formula of vector crossover and mutation, and then keeping whichever candidate solution has the best score or fitness on the optimization problem at hand.



# Bibliography

* Storn, R., Price, K. Differential Evolution – A Simple and Efficient Heuristic for global Optimization over Continuous Spaces. Journal of Global Optimization 11, 341–359 (1997). https://doi.org/10.1023/A:1008202821328 
