#!/usr/bin/env python
# coding: utf-8

import configparser
import logging
import os

from src.utils import get_translation_max

MAIN_PATH = os.getcwd()


class EvodockConfig:
    def __init__(self, ini_file):
        logger = logging.getLogger("evodock.config")
        logger.setLevel(logging.INFO)

        config = self.read_config(ini_file)

        # --- DOCKING PARAMS -----------------------------------+

        self.docking_type_option = config["Docking"].get("type")

        # --- Input Params -----------------------------+
        if config.has_option("inputs", "pose_input"):
            self.pose_input = MAIN_PATH + config["inputs"].get("pose_input")
        else:
            logger.info("input complex not found. Use 'pose_input' parameter")
            exit()

        if config.has_option("inputs", "native_input"):
            self.native_input = MAIN_PATH + config["inputs"].get("native_input")
        else:
            logger.info("native complex not found. Use 'native_input' parameter")
            exit()

        if self.docking_type_option == "Unbound":
            if config.has_option("inputs", "path_ligands"):
                self.path_ligands = MAIN_PATH + config["inputs"].get("path_ligands")
            else:
                logger.info(
                    "path for ligand flexible BB not found. Use 'path_ligands' parameter"
                )
                exit()

            if config.has_option("inputs", "path_receptors"):
                self.path_receptors = MAIN_PATH + config["inputs"].get("path_receptors")
            else:
                logger.info(
                    "path for receptor flexible BB not found. Use 'path_receptors' parameter"
                )
                exit()

            if config.has_option("Docking", "bb_strategy"):
                self.bb_strategy = config["Docking"].get("bb_strategy")
            else:
                logger.info("BB strategy not found. Use 'bb_strategy' parameter")
                exit()

        # --- DE PARAMS -----------------------------------+

        self.scheme = config["DE"].get("scheme")
        self.popsize = config["DE"].getint("popsize")
        self.mutate = config["DE"].getfloat("mutate")
        self.recombination = config["DE"].getfloat("recombination")
        self.maxiter = config["DE"].getint("maxiter")

        # --- MEMETIC PARAMS -----------------------------------+
        if config.has_option("DE", "local_search"):
            self.local_search_option = config["DE"].get("local_search")
        else:
            logger.info("DANGER: local_search is None")
            self.local_search_option = "None"

        # --- OUTPUT FILE -----------------------------------+
        self.jobid = config["outputs"].get("output_file")

        # --- CONFIG STRUCTURE -----------------------------------+
        self.config = config

    def get_max_translation(self):
        # --- Position Params -----------------------------+
        if self.config.has_option("position", "trans_max_magnitude"):
            self.trans_max_magnitude = self.config["position"].getint(
                "trans_max_magnitude"
            )
        else:
            self.trans_max_magnitude = get_translation_max(self.pose_input)
        return self.trans_max_magnitude

    def read_config(self, ini_file):
        config = configparser.ConfigParser()
        config.read(ini_file, encoding="utf-8-sig")
        return config
