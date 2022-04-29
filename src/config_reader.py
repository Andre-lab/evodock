#!/usr/bin/env python
# coding: utf-8

import configparser
import logging
import os

from src.utils import get_translation_max
from itertools import islice

MAIN_PATH = os.getcwd()


class EvodockConfig:
    def __init__(self, ini_file):
        logger = logging.getLogger("evodock.config")
        logger.setLevel(logging.INFO)

        config = self.read_config(ini_file)

        # --- Input Params -----------------------------+
        if config.has_option("inputs", "pose_input"):
            self.pose_input = MAIN_PATH + config["inputs"].get("pose_input")
        else:
            logger.info("input complex not found. Use 'pose_input' parameter")
            exit()

        # --- Symmetrical information ---------------------------+
        self.syminfo = self.extract_syminfo(config)

        # --- DE PARAMS -----------------------------------+
        self.scheme = config["DE"].get("scheme")
        self.popsize = config["DE"].getint("popsize")
        self.mutate = config["DE"].getfloat("mutate")
        self.recombination = config["DE"].getfloat("recombination")
        self.maxiter = config["DE"].getint("maxiter")
        # TODO: DELETE??
        # self.initialization = config["DE"].get("initialization")
        # assert self.initialization in ("gauss", "uniform", None), "initialization can only be of 'gauss' and 'uniform'"

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

        # -- PYMOL ----------------------------------------------+
        if config.has_option("pymol", "on"):
            self.pymol_on = config.getboolean("pymol", "on")
        else:
            self.pymol_on = False
        if config.has_option("pymol", "history"):
            self.pymol_history = config.getboolean("pymol", "history")
        else:
            self.pymol_history = False
        if config.has_option("pymol", "show_local_search"):
            self.show_local_search = config.getboolean("pymol", "show_local_search")
        else:
            self.show_local_search = False

        # -- LOCAL SEARCH OPTIONS -------------------------------+
        if config.has_option("DE", "slide"):
            self.slide = config.getboolean("DE", "slide")
        else:
            self.slide = True

    def extract_syminfo(self, config) -> dict or None:
        """Extracts the symmetrical information from the config file."""
        syminfo = {}
        if config.has_option("inputs", "symmetry_file"):
            syminfo["file"] = MAIN_PATH + config.get("inputs", "symmetry_file")
            if config.has_option("symmetry", "symdofs"):
                syminfo["jumps_str"] = [i.split(":")[0] for i in config.get("symmetry", "symdofs").split(",")]
                syminfo["dofs_str"] = [i.split(":")[1:] for i in config.get("symmetry", "symdofs").split(",")]
                bounds = iter([i for i in config.get("symmetry", "symbounds").split(",")])
                # for nicer mapping we have to map the bounds the same ways as the dofs
                syminfo["bounds"] = [list(islice(bounds, l)) for l in [len(i) for i in syminfo["dofs_str"]]]
                assert len(syminfo["jumps_str"]) == len(syminfo["dofs_str"])
                assert len(syminfo["jumps_str"]) == len(syminfo["bounds"])
            else: # apply defaults
                raise NotImplementedError
        return syminfo if syminfo else None # Return None if empty

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
