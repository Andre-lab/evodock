#!/usr/bin/env python
# coding: utf-8

import logging
import os
import sys

from pyrosetta import init
from pyrosetta.rosetta.core.scoring import CA_rmsd

from src.config_reader import EvodockConfig
from src.differential_evolution import DifferentialEvolutionAlgorithm as DE

from src.options import build_rosetta_flags
from src.population import ScorePopulation
from src.scfxn_fullatom import FAFitnessFunction
from src.single_process import SingleProcessPopulCalculator
from src.utils import get_pose_from_file, get_position_info

MAIN_PATH = os.getcwd()

logging.basicConfig(level=logging.ERROR)


def main():
    config = EvodockConfig(sys.argv[-1])

    pose_input = config.pose_input

    init(extra_options=build_rosetta_flags(config))

    logger = logging.getLogger("evodock")
    logger.setLevel(logging.INFO)

    # --- INIT OPTIONS ------------------------------
    docking_type_option = config.docking_type_option
    jobid = config.out_path

    # --- INIT PROTEIN STRUCTURES -------------------
    native_input = config.native_input
    input_pose = get_pose_from_file(pose_input)
    native_pose = get_pose_from_file(native_input)
    native_pose.conformation().detect_disulfides()
    input_pose.conformation().detect_disulfides()

    # ---- INIT SCORE FUNCTION ------------------------------

    scfxn = FAFitnessFunction(input_pose, native_pose, config)

    # ---- PRINT INIT INFORMATION ---------------------

    position_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(input_pose)]
    )
    native_position_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(native_pose)]
    )

    native_score = scfxn.scfxn_rosetta.score(native_pose)
    input_score = scfxn.scfxn_rosetta.score(input_pose)

    logger.info("==============================")
    logger.info(" input information ")
    logger.info(" input position: " + position_str)
    logger.info(" input pose score {:.2f}".format(input_score))
    logger.info(" native pose score {:.2f}".format(native_score))
    logger.info(" native position: " + native_position_str)
    logger.info(" input vs native rmsd: " + str(CA_rmsd(input_pose, native_pose)))
    logger.info("==============================")

    # ---- START ALGORITHM ---------------------------------
    alg = DE(config, scfxn)
    init_population = alg.init_population(config.popsize, docking_type_option)

    # --- RUN ALGORITHM -------------------------------------

    logger.info("==============================")
    logger.info(" starts EvoDOCK : evolutionary docking process")
    best_pdb = alg.main(init_population)

    # ---- OUTPUT -------------------------------------------
    if config.out_pdb:
        name = jobid + "/final_docked_evo.pdb"
        best_pdb.dump_pdb(name)


if __name__ == "__main__":
    main()
