#!/usr/bin/env python
# coding: utf-8

import logging
import os
import sys

from pyrosetta import init

from src.config_reader import EvodockConfig
from src.differential_evolution import DifferentialEvolutionAlgorithm as DE
from src.differential_evolution import FlexbbDifferentialEvolution as FlexbbDE

from src.options import build_rosetta_flags
from src.scfxn_fullatom import FAFitnessFunction
from src.single_process import SingleProcessPopulCalculator as PopulCalculator
from src.utils import get_pose_from_file, get_position_info, get_symmetric_genotype_str
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from src.utils import get_pose_from_file, get_position_info, get_starting_poses
from pyrosetta.rosetta.core.scoring import CA_rmsd, CA_rmsd_symmetric

MAIN_PATH = os.getcwd()

logging.basicConfig(level=logging.ERROR)


def print_init_information(logger, scfxn, native_pose, input_pose, syminfo: dict = None):
    position_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(input_pose)]
    )
    native_position_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(native_pose)]
    )

    native_score = scfxn.scfxn_rosetta.score(native_pose)
    input_score = scfxn.scfxn_rosetta.score(input_pose)
    input_vs_native_rmsd = CA_rmsd(input_pose, native_pose)

    logger.info("==============================")
    logger.info(" input information ")
    if syminfo:
        logger.info(f" Symmetry has been detected (with genotype: {get_symmetric_genotype_str(native_pose)})")
    logger.info(f" Input position: {position_str}")
    logger.info(f" Input pose score {input_score:.2f}")
    logger.info(f" Native pose score {native_score:.2f}")
    logger.info(f" Native position: {native_position_str}")
    logger.info(f" Input vs native rmsd: {input_vs_native_rmsd:.2f}")
    logger.info("==============================")
# >>>>>>> main


def main():
    config = EvodockConfig(sys.argv[-1])

    init(extra_options=build_rosetta_flags(config))

    logger = logging.getLogger("evodock")
    logger.setLevel(logging.INFO)

    # --- INIT PROTEIN STRUCTURES -------------------
    pose_input = config.pose_input
    native_input = config.native_input
    # input_pose = get_pose_from_file(pose_input)
    # native_pose = get_pose_from_file(native_input)

    input_pose, native_pose = get_starting_poses(pose_input, native_input, config)

    # ---- INIT SCORE FUNCTION ------------------------------
    scfxn = FAFitnessFunction(input_pose, native_pose, config)

    # ---- PRINT INIT INFORMATION ---------------------
    print_init_information(logger, scfxn, native_pose, input_pose, config.syminfo)

    # ---- START ALGORITHM ---------------------------------
    if config.docking_type_option == "Flexbb":
        alg = FlexbbDE(config, scfxn)
    else:
        alg = DE(config, scfxn)

    alg.init_population()

    # --- RUN ALGORITHM -------------------------------------
    logger.info("==============================")
    logger.info(" starts EvoDOCK : evolutionary docking process")
    best_pdb = alg.main()

    # ---- OUTPUT -------------------------------------------
    logger.info(" end EvoDOCK")
    logger.info("==============================")
    if config.out_pdb:
        name = config.out_path + "/final_docked_evo.pdb"
        best_pdb.dump_pdb(name)

if __name__ == "__main__":
    main()
