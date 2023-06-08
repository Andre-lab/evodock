#!/usr/bin/env python
# coding: utf-8


class GreedySelection:
    def apply(self, trials, population, selection):
        if selection == "total":
            gen_scores = [target.score for target in population]
            trial_scores = [ind.score for ind in trials]
        else: #selection ==  "interface"
            gen_scores = [target.i_sc for target in population]
            trial_scores = [ind.i_sc for ind in trials]

        for j in range(0, len(population)):
            target_idx = j
            score_trial = trial_scores[j]
            if selection == "total":
                score_target = population[target_idx].score
            else: #selection ==  "interface"
                score_target = population[target_idx].i_sc
            if score_trial < score_target:
                population[target_idx] = trials[j]
                gen_scores[target_idx] = score_trial
            else:
                gen_scores[target_idx] = score_target

        return population, gen_scores, trial_scores
