#!/usr/bin/env python
# coding: utf-8


import logging
import os
import sys

from pyrosetta import init

from src.config_reader import EvodockConfig
from src.differential_evolution import DifferentialEvolutionAlgorithm as DE
from src.options import init_global_docking, init_local_docking
from src.population import ScorePopulation
from src.scfxn_fullatom import FAFitnessFunction
from src.single_process import SingleProcessPopulCalculator as PopulCalculator
from src.utils import get_pose_from_file, get_position_info

MAIN_PATH = os.getcwd()

logging.basicConfig(level=logging.ERROR)


def main():
    config = EvodockConfig(sys.argv[-1])

    pose_input = config.pose_input
    if config.docking_type_option == "Local":
        init(extra_options=init_local_docking(pose_input))
    else:
        init(extra_options=init_global_docking(pose_input))

    logger = logging.getLogger("evodock")
    logger.setLevel(logging.INFO)

    # --- Position Params -----------------------------+
    trans_max_magnitude = config.get_max_translation()

    # --- DOCKING PARAMS -----------------------------------+
    docking_type_option = config.docking_type_option

    # --- LS PARAMS -----------------------------------+
    local_search_option = config.local_search_option

    # --- OUTPUT --------------------------------------+
    jobid = config.jobid
    if "/" in jobid:
        os.makedirs("/".join(jobid.split("/")[:-1]), exist_ok=True)

    # --- INIT -----------------------------------------+
    native_pose = get_pose_from_file(pose_input)
    scfxn = FAFitnessFunction(native_pose, trans_max_magnitude)

    position_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(native_pose)]
    )
    native_score = scfxn.scfxn_rosetta.score(native_pose)

    score_popul = ScorePopulation(scfxn, jobid, local_search_option)
    popul_calculator = PopulCalculator(score_popul)

    logger.info("==============================")
    logger.info(" native information ")
    logger.info(" native position: " + position_str)
    logger.info(" native pose score {:.2f}".format(native_score))

    logger.info("==============================")
    alg = DE(popul_calculator, config)
    init_population = alg.init_population(config.popsize, docking_type_option)

    # --- RUN -----------------------------------------+
    logger.info("==============================")
    logger.info(" starts EvoDOCK : evolutionary docking process")
    population = popul_calculator.run(init_population)
    _, best_pdb = alg.main(population)
    popul_calculator.terminate()
    name = jobid.replace(".log", "_final_docked_evo.pdb")
    best_pdb.dump_pdb(name)


if __name__ == "__main__":
    main()
