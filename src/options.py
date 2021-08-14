#!/usr/bin/env python
# coding: utf-8


def init_global_docking(filename):
    opts = [
        "-mute all",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-include_current True",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-unboundrot {}".format(filename),
    ]
    return " ".join(opts)


def init_local_docking(filename):
    opts = [
        "-mute all",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-dock_pert 3 8",
        "-docking:docking_centroid_inner_cycles 1",
        "-docking:docking_centroid_outer_cycles 1",
        "-include_current True",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-unboundrot {}".format(filename),
    ]
    return " ".join(opts)