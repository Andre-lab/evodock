#!/usr/bin/env python
# coding: utf-8
import random
import numpy as np


def build_rosetta_flags(config):
    filename = config.pose_input
    unboundrot = filename
    # create unbound rotamers from the template
    if config.template is not None:
        unboundrot = config.template
    # create unbound rotamers from the
    elif config.syminfo is None and config.flexbb:
        unboundrot = f"{filename[0]} {filename[1]}" if isinstance(filename, list) else f"{filename}",
    seed = config.seed
    # add fixed options
    opts = [
        "-mute all",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-include_current True",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-extrachi_cutoff 1",
        "-unboundrot " + unboundrot,
    ]
    # add additional_rosetta_options
    if config.rosetta_options is not None:
        for k,v in config.rosetta_options:
            opts.append(f"-{k} {v}")
    # add seed if set
    if seed:
        if len(str(seed)) > 10:
            exit("Seed is to High. Use 10 or less digits.")
        opts += ["-run:constant_seed", f"-run:jran {seed}"]
        random.seed(seed)
        np.random.seed(seed)
    return " ".join(opts)


# def init_local_docking(filename, syminfo: dict = None):
#     raise NotImplementedError("Shouldnt use this. Need to be refactored!!")
#     opts = [
#         "-mute all",
#         # "-out:level 9999",
#         "-output_virtual True",
#         "-detect_disulf True",
#         "-docking:dock_mcm_first_cycles 1",
#         "-docking:dock_mcm_second_cycles 1",
#         "-dock_pert 3 8",
#         "-relax:jump_move true",
#         "-relax:bb_move true",
#         "-relax:chi_move true",
#         "-docking:docking_centroid_inner_cycles 1",
#         "-docking:docking_centroid_outer_cycles 1",
#         "-include_current True",
#         "-ex1",
#         "-ex2aro",
#         "-extrachi_cutoff 1",
#         "-use_input_sc",
#         "-unboundrot {}".format(filename),
#     ]
#     return " ".join(opts) + add_syminfo_to_init(syminfo)
