#!/usr/bin/env python
# coding: utf-8


import logging
import os
import sys

from pyrosetta import init

from src.config_reader import EvodockConfig
from src.differential_evolution import DifferentialEvolutionAlgorithm as DE
from src.population import ScorePopulation
from src.scfxn_fullatom import FAFitnessFunction
from src.single_process import SingleProcessPopulCalculator as PopulCalculator
from src.utils import get_pose_from_file, get_position_info, get_symmetric_genotype_str
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric

MAIN_PATH = os.getcwd()

logging.basicConfig(level=logging.ERROR)


def init_options_fa(filename):
    opts = [
        "-mute all",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-include_current True",
        "-initialize_rigid_body_dofs 1",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-unboundrot {}".format(filename),
    ]
    return " ".join(opts)


def main():
    config = EvodockConfig(sys.argv[-1])

    pose_input = config.pose_input

    init(extra_options=init_options_fa(pose_input))

    logger = logging.getLogger("evodock")
    logger.setLevel(logging.INFO)

    # --- Position Params -----------------------------+
    trans_max_magnitude = config.get_max_translation()

    # --- LS PARAMS -----------------------------------+
    local_search_option = config.local_search_option

    # --- OUTPUT --------------------------------------+
    jobid = config.jobid

    # --- Symmetry information ------------------------+
    syminfo = config.syminfo

    # --- INIT ----------------------------------------+
    native_pose = get_pose_from_file(pose_input, syminfo)
    scfxn = FAFitnessFunction(native_pose, trans_max_magnitude, syminfo)

    position_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(native_pose)]
    )
    native_score = scfxn.scfxn_rosetta.score(native_pose)

    score_popul = ScorePopulation(scfxn, jobid, local_search_option)
    popul_calculator = PopulCalculator(score_popul, syminfo)

    if is_symmetric(native_pose):
        logger.info("==============================")
        logger.info(" Running EvoDOCK with symmetry")
        logger.info(f" The genotype is: {get_symmetric_genotype_str(native_pose)}")

    logger.info("==============================")
    logger.info(" native information ")
    logger.info(" native position: " + position_str)
    logger.info(" native pose score {:.2f}".format(native_score))

    logger.info("==============================")
    alg = DE(popul_calculator, config)
    init_population = alg.init_population()

    # --- RUN -----------------------------------------+
    logger.info("==============================")
    logger.info(" starts EvoDOCK : evolutionary docking process")
    population = popul_calculator.run(init_population)
    alg.main(population)
    popul_calculator.terminate()

if __name__ == "__main__":
    main()
