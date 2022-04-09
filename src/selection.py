#!/usr/bin/env python
# coding: utf-8


class GreedySelection:
    def apply(self, trials, population):
        gen_scores = [target.score for target in population]
        trial_scores = [ind.score for ind in trials]

        for j in range(0, len(population)):
            target_idx = j
            score_trial = trial_scores[j]
            score_target = population[target_idx].score
            if score_trial < score_target:
                population[target_idx] = trials[j]
                gen_scores[target_idx] = score_trial
            else:
                gen_scores[target_idx] = score_target

        return population, gen_scores, trial_scores
