#!/usr/bin/env python
# coding: utf-8

import random

from src.scfxn_fullatom import l2_norm


def nearest_neighbor(population, trial, near_option="TRIAL"):
    distance_values = []
    for i in range(0, len(population)):
        if near_option == "TRIAL":
            if random.uniform(0, 1) < 0.1:
                dst = l2_norm(population[i].genotype[3:], trial[3:])
                distance_values.append((i, dst))
        else:
            dst = l2_norm(population[i].genotype, trial)
            distance_values.append((i, dst))
    distance_values.sort(key=lambda x: x[1])
    if near_option == "TRIAL":
        if len(distance_values) < 1:
            nearest = int(random.uniform(0, len(population)))
        else:
            nearest = distance_values[0][0]
    else:
        prune = [x for x in distance_values if x[1] != 0]
        nearest = prune[0][0]
    return nearest


class GreedySelection:
    def apply(self, trials, population):
        gen_scores = [target.score for target in population]
        trial_scores = [ind.score for ind in trials]
        for j in range(0, len(population)):
            # target_idx = self.nearest_neighbor(population, trials[j].genotype)
            target_idx = j
            score_trial = trial_scores[j]
            score_target = population[target_idx].score

            if score_trial < score_target:
                population[target_idx] = trials[j]
                gen_scores[target_idx] = score_trial
            else:
                gen_scores[target_idx] = score_target

        return population, gen_scores, trial_scores


class CrowdingSelection:
    def apply(self, trials, population):
        gen_scores = [target.score for target in population]
        trial_scores = [ind.score for ind in trials]
        for j in range(0, len(population)):
            target_idx = nearest_neighbor(population, trials[j].genotype)
            score_trial = trial_scores[j]
            score_target = population[target_idx].score

            if score_trial < score_target:
                population[target_idx] = trials[j]
                gen_scores[target_idx] = score_trial
            else:
                gen_scores[target_idx] = score_target

        return population, gen_scores, trial_scores


class EliteSelection:
    def __init__(self, popul_calculator):
        self.popul_calculator = popul_calculator
        self.ind_size = popul_calculator.cost_func.size()
        self.bounds = [(-1, 1)] * self.ind_size

    def apply(self, trials, population):
        score_popul = [(p.score, p) for p in population]
        score_trials = [(p.score, p) for p in trials]
        total_popul_gen = score_popul + score_trials
        total_popul_gen.sort(key=lambda x: x[0], reverse=False)

        elite_size = int(len(population) * 0.25)
        top_individuals = total_popul_gen[:elite_size]
        top_individuals = [x[1] for x in top_individuals]

        restart_popul = []
        b_gen = top_individuals[0].genotype
        #  => 60 /180
        bound_range = 0.35
        bounds_rotation = [
            (b_gen[j] - bound_range, b_gen[j] + bound_range) for j in range(0, 3)
        ]
        # => 5 /180
        bound_range = 0.03
        bounds_translation = [
            (b_gen[j] - bound_range, b_gen[j] + bound_range)
            for j in range(3, len(b_gen))
        ]
        bounds = bounds_rotation + bounds_translation
        for i in range(0, len(population) - len(top_individuals)):
            indv = []
            for j in range(len(self.bounds)):
                indv.append(random.uniform(bounds[j][0], bounds[j][1]))
            restart_popul.append(indv)

        restart_popul[-1] = [
            -126.20 / 180,
            59.54 / 180,
            -43.19 / 180,
            30.59 / 70,
            -33.25 / 70,
            -29.16 / 70,
        ]

        score_restarted_inds = self.popul_calculator.run(restart_popul)

        next_gen_popul = top_individuals + score_restarted_inds
        gen_scores = [target.score for target in next_gen_popul]
        trial_scores = [ind.score for ind in trials]
        return next_gen_popul, gen_scores, trial_scores
