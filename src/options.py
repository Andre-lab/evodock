#!/usr/bin/env python
# coding: utf-8


def build_rosetta_flags(config):
    if config.docking_type_option in ["Local", "Flexbb", "Refine"]:
        return init_local_docking(config.pose_input, config.syminfo)
    else:
        return init_global_docking(config.pose_input, config.syminfo)

def add_syminfo_to_init(syminfo: dict = None):
    #         "-unmute protocols.simple_moves_symmetry.SymDockingInitialPerturbation",
    #         "-out:file:output_pose_energies_table false", # FIXME bypasses this error: Energies::residue_total_energies( int const seqpos ): variable seqpos is out of range!
    if syminfo and syminfo.get("initialize_rigid_body_dofs"):
        return " -initialize_rigid_body_dofs 1 "
    else:
        return ""

def init_global_docking(filename, syminfo: dict = None):
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
    return " ".join(opts) + add_syminfo_to_init(syminfo)


def init_local_docking(filename, syminfo: dict = None):
    opts = [
        "-mute all",
        # "-out:level 9999",
        "-output_virtual True",
        "-detect_disulf True",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-dock_pert 3 8",
        "-relax:jump_move true",
        "-relax:bb_move true",
        "-relax:chi_move true",
        "-docking:docking_centroid_inner_cycles 1",
        "-docking:docking_centroid_outer_cycles 1",
        "-include_current True",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-unboundrot {}".format(filename),
    ]
    return " ".join(opts) + add_syminfo_to_init(syminfo)
