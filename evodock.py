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
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from cubicsym.cubicsetup import CubicSetup
from cubicsym.assembly.cubicassembly import CubicSymmetricAssembly
from pathlib import Path

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

def output_models(alg, config, logger):
    """Output the models. Does KMeans clustering if set to True."""
    vals = alg.all_docks["genotype"]
    n_models = config.n_models
    # check that the individuals are all different
    # assert not sum([v is vv for v in vals for vv in vals]) - 1 < len(vals)
    # Cluster or not and make into a dataframe
    df = pd.DataFrame(alg.all_docks)
    if config.cluster:
        # do kmeans clustering
        kmeans = KMeans(n_clusters=n_models).fit(vals)
        df["cluster"] = kmeans.labels_
        if config.selection == "total":
            df = df.sort_values("score").groupby("cluster").first()
        else:
            df = df.sort_values("i_sc").groupby("cluster").first()
    else:
        if config.selection == "total":
            df = df.sort_values("score")[0:n_models]
        else:
            df = df.sort_values("i_sc")[0:n_models]
    # output the models
    logger.info("==============================")
    # make a subfolder
    output_folder = config.out_path + "/structures"
    Path(output_folder).mkdir(exist_ok=True)
    for n, ind in enumerate(df["ind"].values, 1):
        if config.flexbb:
            pose, _, _, _ = alg.scfxn.local_search.local_search_strategy.get_pose(ind)
        else:
            pose = alg.scfxn.local_search.local_search_strategy.get_pose(ind)
        name = output_folder + f"/final_evodock_unrelaxed_ranked_{n}{{}}"
        if is_symmetric(pose):
            logger.info(" Outputting the following models:")
            cs = CubicSetup(config.syminfo.input_symdef)
            if cs.is_cubic():
                # output full symmetric structure
                name_full = name.format(f"_full.cif")
                logger.info(f" {name_full}")
                CubicSymmetricAssembly.from_pose_input(pose, cs).output(filename=name_full, format="cif")
                # output symmetry file
                name_symm = name.format(f".symm")
                logger.info(f" {name_symm}")
                cs.output(name_symm)
                # output input pdb
                name_input = name.format(f"_INPUT.pdb")
                logger.info(f" {name_input}")
                cs.make_asymmetric_pose(pose).dump_pdb(name_input)
            else:
                raise NotImplementedError("Only Cubic symmetry allowed")
    logger.info("==============================")

def main():
    #========
    print("="*43)
    print(" "*18, "EvoDOCK", " "*18)
    print("="*43)

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
    logger.info(" Starts EvoDOCK")
    alg.main()

    # ---- OUTPUT -------------------------------------------
    output_models(alg, config, logger)

    logger.info(" End EvoDOCK")
    logger.info("==============================")


if __name__ == "__main__":
    main()
