#!/usr/bin/env python
# coding: utf-8

import configparser
import logging
import os
import sys

import numpy as np
# from mpi4py import MPI
from pyrosetta import SwitchResidueTypeSetMover, init
from pyrosetta.rosetta.protocols.simple_moves import ReturnSidechainMover

from src.differential_evolution import DifferentialEvolutionAlgorithm as DE
from src.init_random_positions import start_input_poses
from src.population import ScorePopulation
from src.scfxn_fullatom import FAFitnessFunction
from src.scfxn_lowres import LowResFitnessFunction
# from src.mpi_utils import MasterProcess, Worker
from src.single_process import SingleMasterProcess as MasterProcess
# from src.scfxn_zernike import FitnessFunction
from src.utils import get_position_info

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()


MAIN_PATH = os.getcwd()

logging.basicConfig(level=logging.ERROR)
# logging.disable(logging.INFO)


def init_options(reference_input):
    opts = [
        "-mute all",
        "-docking:randomize1",
        "-docking:randomize2",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-docking:dock_mcm_trans_magnitude 2",
        "-docking:dock_mcm_rot_magnitude 3",
        "-docking:low_res_protocol_only True",
        "-docking:no_filters",
        "-zernike_descriptor:zernike_descriptor_file {}".format(reference_input),
        "-zernike_descriptor::grid_size 64",
        "-docking::docklowres_trans_magnitude 3",
        "-docking::docklowres_rot_magnitude 5",
        "-ex1",
        "-ex2aro",
    ]
    return " ".join(opts)


def init_options_fa(filename):
    opts = [
        "-mute all",
        # "-mute all -unmute core.pack.rotamers core.pack.rotamer_trials core.pack.task protocols protocols.docking.DockMCMProtocol protocols.docking.DockMCMCycle protocols.moves.MoverContainer",
        # "-out:level 1000",
        # "-docking:dock_mcm_trans_magnitude 0",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-dock_pert 3 8",
        "-docking:docking_centroid_inner_cycles 1",
        "-docking:docking_centroid_outer_cycles 1",
        # "-docking_local_refine",
        "-include_current True",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        # "-extrachi_cutoff 1",
        "-unboundrot {}".format(filename),
        # "-mh:score:use_ss1 false",
        # "-mh:score:use_ss2 false",
        # "-mh:score:use_aa1 true",
        # "-mh:score:use_aa2 true",
    ]
    return " ".join(opts)


def read_config(ini_file):
    config = configparser.ConfigParser()
    config.read(ini_file, encoding="utf-8-sig")
    return config


def get_translation_max(native_pose, dock_pose):
    jump_num = 1
    flexible_jump = dock_pose.jump(jump_num)
    translation = np.asarray(flexible_jump.get_translation())
    return max(translation) + 10


def main():
    config = read_config(sys.argv[-1])
    if config.has_option("inputs", "reference_input"):
        reference_input = MAIN_PATH + config["inputs"].get("reference_input")
    else:
        reference_input = ""

    pose_input = MAIN_PATH + config["inputs"].get("pose_input")

    # init(extra_options=init_options(reference_input))
    init(extra_options=init_options_fa(pose_input))

    logger = logging.getLogger("evodock")
    logger.setLevel(logging.INFO)

    # --- Position Params -----------------------------+
    trans_max_magnitude = config["position"].getint("trans_max_magnitude")

    # --- DE PARAMS -----------------------------------+

    scheme = config["DE"].get("scheme")
    popsize = config["DE"].getint("popsize")
    mutate = config["DE"].getfloat("mutate")
    recombination = config["DE"].getfloat("recombination")
    maxiter = config["DE"].getint("maxiter")

    if config.has_option("DE", "local_search"):
        local_search_option = config["DE"].get("local_search")
    else:
        logger.info("DANGER: local_search is None")
        local_search_option = "None"

    if config.has_option("protocol", "refinement"):
        refinement_option = config["protocol"].getboolean("refinement")
    else:
        logger.info("DANGER: refinement is False by default")
        # kind of deprecated: perhaps for a next version
        refinement_option = False

    if config.has_option("protocol", "include_lowres_stage"):
        include_lowres_stage_option = config["protocol"].getboolean(
            "include_lowres_stage"
        )
    else:
        logger.info("DANGER: include_lowres_stage is False by default")
        # kind of deprecated: perhaps for a next version
        include_lowres_stage_option = False

    # --- INIT -----------------------------------------+
    jobid = config["outputs"].get("output_file")
    native_pose, init_state_pose = start_input_poses(pose_input)
    # -- TESTING PURPOUSES ----
    position_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(native_pose)]
    )
    input_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(init_state_pose)]
    )

    # -- TESTING PURPOUSES ----

    refinement_option = True
    trans_max_magnitude = get_translation_max(native_pose, init_state_pose)
    scfxn = FAFitnessFunction(
        native_pose, init_state_pose, trans_max_magnitude, refinement_option
    )

    score_popul = ScorePopulation(scfxn, jobid, local_search_option)
    rank = 0
    size = 1
    master_calculator = MasterProcess(size, score_popul)

    if rank == 0:
        logger.info("==============================")
        logger.info(" init the params ")
        logger.info(" native_position: " + position_str)
        logger.info(" input_position: " + input_str)
        alg = DE(
            master_calculator, scheme, popsize, mutate, recombination, maxiter, jobid,
        )
        init_population = alg.init_population()
        # --- RUN -----------------------------------------+
        refinement_option = False
        scfxn = FAFitnessFunction(
            native_pose, init_state_pose, trans_max_magnitude, refinement_option
        )

        score_popul = ScorePopulation(scfxn, jobid, local_search_option)
        master_calculator = MasterProcess(size, score_popul)
        alg = DE(
            master_calculator, scheme, popsize, mutate, recombination, maxiter, jobid,
        )

        logger.info("==============================")
        logger.info(" starts high res evolution")
        logger.info(" native pose {}".format(scfxn.scfxn_rosetta.score(native_pose)))
        high_res_population = master_calculator.run(init_population)
        alg.main(high_res_population)
        master_calculator.terminate()
    # else:
    #     Worker(rank, size, score_popul).run()


if __name__ == "__main__":
    main()
