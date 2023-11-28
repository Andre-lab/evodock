#!/usr/bin/env python
# coding: utf-8

import logging
import os
import sys

import pandas as pd
from pyrosetta import init
from src.config_reader import EvodockConfig
from src.differential_evolution import DifferentialEvolutionAlgorithm as DE
from src.differential_evolution import FlexbbDifferentialEvolution as FlexbbDE
from src.options import build_rosetta_flags
from src.scfxn_fullatom import FAFitnessFunction
from src.utils import get_position_info, initialize_starting_poses
from src.dock_metric import DockMetric, CubicDockMetric
from sklearn.cluster import KMeans

MAIN_PATH = os.getcwd()

logging.basicConfig(level=logging.ERROR)


def print_init_information(logger, scfxn, native_pose, input_pose, dockmetric, syminfo: dict = None, native_symmetric_pose=None):
    input_score = f"{scfxn.scfxn_rosetta.score(input_pose) :.2f}"
    position_str = ", ".join(
        ["{:.2f}".format(e) for e in get_position_info(input_pose, syminfo)]
    )
    if native_pose is not None:
        native_score = f"{scfxn.scfxn_rosetta.score(native_symmetric_pose if native_symmetric_pose is not None else native_pose) :.2f}"
        native_position_str = ", ".join(
            ["{:.2f}".format(e) for e in get_position_info(native_symmetric_pose if native_symmetric_pose is not None else native_pose, syminfo)]
        )
        input_vs_native_rmsd = dockmetric.ca_rmsd(input_pose)

    logger.info("==============================")
    logger.info(" Input information ")
    logger.info(f" Input position: {position_str}{' (is symmetrical)' if syminfo else ''}")
    logger.info(f" Input pose score {input_score}")
    if native_pose is not None:
        logger.info(f" Native pose score {native_score}")
        logger.info(f" Native position: {native_position_str}")
        logger.info(f" Input vs native rmsd: {input_vs_native_rmsd:.2f}")
    logger.info("==============================")

def initialize_dock_metric(config, native, input_pose):
    if config.syminfo:
        ### fixme: delete this
        jump_ids, dof_ids, trans_mags = [], [], []
        for a, b, c in config.syminfo.normalize_trans_map:
            jump_ids.append(a)
            dof_ids.append(b)
            trans_mags.append(c)
        ###
        return CubicDockMetric(native, input_pose, config.syminfo.native_symdef, config.syminfo.input_symdef,
                               jump_ids=jump_ids, dof_ids=dof_ids, trans_mags=trans_mags, use_map=config.rmsd_map)
    else:
        return DockMetric(native)

def output_best_clusters(alg, config, logger):
    # do kmeans clustering and output the final results.
    vals = alg.all_docks["genotype"]
    kmeans = KMeans(n_clusters=5, random_state=0).fit(vals)
    alg.all_docks["cluster"] = kmeans.labels_
    df = pd.DataFrame(alg.all_docks)
    if config.selection == "total":
        df = df.sort_values("score").groupby("cluster").first()
    else:
        df = df.sort_values("i_sc").groupby("cluster").first()
    logger.info("==============================")
    logger.info(" Clustering and outputting the following models:")
    for n, ind in enumerate(df["ind"].values, 1):
        if config.flexbb:
            pose, _, _, _ = alg.scfxn.local_search.local_search_strategy.get_pose(ind)
        else:
            pose = alg.scfxn.local_search.local_search_strategy.get_pose(ind)
        name = config.out_path + f"/final_evodock_unrelaxed_ranked_{n}.pdb"
        logger.info(f" {name}")
        pose.dump_pdb(name)
    logger.info("==============================")

def main():
    config = EvodockConfig(sys.argv[-1])

    init(extra_options=build_rosetta_flags(config))

    logger = logging.getLogger("evodock")
    logger.setLevel(logging.INFO)

    # --- INIT PROTEIN STRUCTURES -------------------
    input_pose, native_pose, native_symmetric_pose = initialize_starting_poses(config)

    # --- INIT METRIC CALCULATOR --------------------
    dockmetric = initialize_dock_metric(config, native_pose, input_pose)

    # ---- INIT SCORE FUNCTION ----------------------
    scfxn = FAFitnessFunction(input_pose, native_pose, config, dockmetric, native_symmetric_pose)

    # ---- PRINT INIT INFORMATION -------------------
    print_init_information(logger, scfxn, native_pose, input_pose, dockmetric, config.syminfo, native_symmetric_pose)

    # ---- START ALGORITHM --------------------------
    if config.flexbb:
        alg = FlexbbDE(config, scfxn)
    else:
        alg = DE(config, scfxn)

    alg.init_population()

    # --- RUN ALGORITHM -------------------------------------
    logger.info("==============================")
    logger.info(" starts EvoDOCK : evolutionary docking process")
    alg.main()


    # ---- OUTPUT -------------------------------------------
    output_best_clusters(alg, config, logger)

    logger.info(" end EvoDOCK")
    logger.info("==============================")



if __name__ == "__main__":
    main()
