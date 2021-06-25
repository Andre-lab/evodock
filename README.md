
# Protein-Protein docking using a Memetic Algorithm: EvoDOCK

:warning: under construction

# Dependencies

* PyRosetta==4
* numpy==1.19.4
* vector3d



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

# Differential Evolution Algorithm

Differential Evolution [Price97] is a population-based search method. DE creates new candidate solutions by combining existing ones according to a simple formula of vector crossover and mutation, and then keeping whichever candidate solution has the best score or fitness on the optimization problem at hand.



# Bibliography

* Storn, R., Price, K. Differential Evolution – A Simple and Efficient Heuristic for global Optimization over Continuous Spaces. Journal of Global Optimization 11, 341–359 (1997). https://doi.org/10.1023/A:1008202821328 
