#!/usr/bin/env python
# coding: utf-8

import random

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pyrosetta import Pose, init, pose_from_file

from src.differential_evolution import Individual
from src.population import ScorePopulation
from src.scfxn_fullatom import FAFitnessFunction, l2_norm
from src.single_process import SingleMasterProcess as MasterProcess


def select_seed(num_population):
    # -- sort popul to select the first seed

    sort_list = sorted(
        list(zip(range(len(num_population)), num_population)),
        key=lambda x: x[1][1].score,
    )
    # num_population = sorted(num_population, key=lambda x: x[1].score)
    for i in range(0, len(sort_list)):
        if sort_list[i][1][0] == -1:
            seed_idx = sort_list[i][0]
            return seed_idx


def calculate_distances(num_population, seed_idx):
    # -- calculate the distance of individuals that were not assigned to any subpopul
    distances_values = [
        (
            i,
            l2_norm(
                num_population[i][1].genotype, num_population[seed_idx][1].genotype
            ),
        )
        for i in range(len(num_population))
        if num_population[i][0] == -1
    ]
    return distances_values


def restart_individual():
    bounds = [(-1, 1)] * 6
    indv = []
    for j in range(len(bounds)):
        indv.append(random.uniform(bounds[j][0], bounds[j][1]))
    return Individual(indv, 1000, 1000)


def assign_subpopul(num_population, count_subpopuls):
    tmp_subpopul = [
        (i, num_population[i][1])
        for i in range(len(num_population))
        if num_population[i][0] == count_subpopuls
    ]

    original_index = [i for i, _ in tmp_subpopul]
    tmp_subpopul = [ind for _, ind in tmp_subpopul]
    index = min(range(len(tmp_subpopul)), key=lambda i: tmp_subpopul[i].score)

    for i in range(len(tmp_subpopul)):
        if (tmp_subpopul[i].score - tmp_subpopul[index].score < 3) and i != index:
            tmp_subpopul[i] = restart_individual()

    return original_index, tmp_subpopul


def create_subpopul(population, size_subpopul, goal_subpopuls):
    count_subpopuls = 0
    subpopulation = {}

    # -- assign to each individual a subpopulation number
    num_population = list(zip([-1] * len(population), population))
    while len(subpopulation) < goal_subpopuls:
        if len(subpopulation) < (goal_subpopuls - 1):
            # -- sort popul to select the first seed
            seed_idx = select_seed(num_population)
            num_population[seed_idx][1] == count_subpopuls
            # -- calculate the distance of individuals that were not assigned to any subpopul
            distances_values = calculate_distances(num_population, seed_idx)

            # -- nearest n elemenets are assigned to the current seed subpopul
            distances_values.sort(key=lambda x: x[1])
            for i in range(size_subpopul):
                num_population[distances_values[i][0]] = (
                    count_subpopuls,
                    num_population[distances_values[i][0]][1],
                )

            orig_idx, subpopulation[count_subpopuls] = assign_subpopul(
                num_population, count_subpopuls
            )
        else:
            orig_idx, subpopulation[count_subpopuls] = assign_subpopul(
                num_population, -1
            )
        count_subpopuls += 1

    return subpopulation


def main():
    init(extra_options="-mute all")
    pose = Pose()
    pose_from_file(pose, "../easy_dock/2hrk_AB.prepack.pdb")

    scfxn = FAFitnessFunction(pose, pose, 50)
    score_popul = ScorePopulation(scfxn, "test", "custom_packer")

    bounds = [(-1, 1)] * 6
    population = []

    master_calculator = MasterProcess(1, score_popul)
    for i in range(0, 1000):
        indv = []
        for j in range(6):
            indv.append(random.uniform(bounds[j][0], bounds[j][1]))
        population.append(Individual(indv, 0, 1000))

    population = master_calculator.run(population)
    goal_subpopuls = 5
    size_subpopul = int(len(population) / goal_subpopuls)
    print(len(population))
    print("subpopulation to be created {} {} ".format(size_subpopul, goal_subpopuls))
    subpopulation = create_subpopul(population, size_subpopul, goal_subpopuls)

    popul_data = []

    for i, vec in subpopulation.items():
        print("subpopulation {}".format(i))
        needs_recalculate = len([k for k in vec if k.score == 1000]) > 0
        if needs_recalculate:
            print("popul needs recalculate")
            subpopulation[i] = master_calculator.run(subpopulation[i])
        for j, ind in enumerate(subpopulation[i]):
            print("ind {} : {} ".format(j, ind.score))
            popul_data.append({"rmsd": ind.rmsd, "energy": ind.score, "subpopul": i})

    df = pd.DataFrame(popul_data)
    ax = sns.scatterplot(
        x="rmsd", y="energy", data=df, hue="subpopul", palette="tab10",
    )
    fig = ax.get_figure()
    fig.set_size_inches(8.7, 5.0)
    fig.tight_layout()
    name = "example_of_subpopuls.png"
    fig.savefig(name)
    plt.close()


if __name__ == "__main__":
    main()
