#!/usr/bin/env python
# coding: utf-8
import random
import numpy as np


def build_rosetta_flags(config):
    # if config.docking_type_option in ["Local", "Flexbb", "Refine"]:
    #     raise NotImplementedError("Shouldnt use this. Need to be refactored!!")
    #     return init_local_docking(config.pose_input, config.syminfo, config.seed)
    # else:
    return init_global_docking(config.pose_input, config.syminfo, config.seed)

def add_syminfo_to_init(syminfo: dict = None):
    #         "-unmute protocols.simple_moves_symmetry.SymDockingInitialPerturbation",
    #         "-out:file:output_pose_energies_table false", # FIXME bypasses this error: Energies::residue_total_energies( int const seqpos ): variable seqpos is out of range!
    if syminfo and syminfo.get("initialize_rigid_body_dofs"):
        return " -initialize_rigid_body_dofs 1 "
    else:
        return ""

def init_global_docking(filename, syminfo: dict = None, seed:int=None):
    opts = [
        "-mute all",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-include_current True",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-extrachi_cutoff 1",
        "-unboundrot {}".format(filename),
    ]
    if seed:
        if len(str(seed)) > 10:
            exit("Seed is to High. Use 10 or less digits.")
        opts += ["-run:constant_seed", f"-run:jran {seed}"]
        random.seed(seed)
        np.random.seed(seed)
    return " ".join(opts) + add_syminfo_to_init(syminfo)


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
