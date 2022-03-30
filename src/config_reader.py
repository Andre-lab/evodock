#!/usr/bin/env python
# coding: utf-8

import configparser
import logging
import os
import ast

from src.utils import get_translation_max

MAIN_PATH = os.getcwd()


class EvodockConfig:
    def __init__(self, ini_file):
        self.p = os.getcwd()
        self.logger = logging.getLogger("evodock.config")
        self.logger.setLevel(logging.INFO)

        config = self.read_config(ini_file)

        self.check_required_parameters(config)
        self.load_required_parameters(config)
        self.load_output_parameters(config)

        # depend on docking type
        self.load_flexible_docking_parameters(config)
        self.load_refine_docking_parameters(config)

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

    def check_required_parameters(self, config):
        # --- REQUIRED PARAMS -----------------------------------+
        req_parameters = [
            ("Docking", "type"),
            ("Inputs", "pose_input"),
            ("Inputs", "native_input"),
            ("DE", "scheme"),
            ("DE", "popsize"),
            ("DE", "mutate"),
            ("DE", "recombination"),
            ("DE", "maxiter"),
            ("DE", "local_search"),
        ]
        all_required_found = True
        for param in req_parameters:
            if config.has_option(param[0], param[1]):
                self.docking_type_option = config[param[0]].get(param[1])
            else:
                self.logger.info(f"[{param[0]}]{param[1]} not found.")
                all_required_found = False
        if not all_required_found:
            exit()

    def load_required_parameters(self, config):
        self.docking_type_option = config["Docking"].get("type")

        # -- INPUT PARAMETERS ---------------------------------------+
        self.pose_input = self.p + config["Inputs"].get("pose_input")
        if not os.path.isfile(self.pose_input):
            self.logger.info(f"input file not found: {self.pose_input}")
            exit()

        self.native_input = self.p + config["Inputs"].get("native_input")
        if not os.path.isfile(self.native_input):
            self.logger.info(f"input file not found: {self.native_input}")
            exit()

        # --- DE PARAMS ---------------------------------------------+
        self.scheme = config["DE"].get("scheme")
        self.popsize = config["DE"].getint("popsize")
        self.mutate = config["DE"].getfloat("mutate")
        self.recombination = config["DE"].getfloat("recombination")
        self.maxiter = config["DE"].getint("maxiter")

        # --- MEMETIC PARAMS -----------------------------------+
        if config.has_option("DE", "local_search"):
            self.local_search_option = config["DE"].get("local_search")
        else:
            self.logger.info("DANGER: local_search is None")
            self.local_search_option = "None"

    def load_flexible_docking_parameters(self, config):
        if self.docking_type_option == "Flexbb":
            if config.has_option("Flexbb", "path_ligands"):
                self.path_ligands = self.p + config["Flexbb"].get("path_ligands")
            else:
                self.logger.info(
                    "path for ligand flexible BB not found. Use 'path_ligands' parameter"
                )
                exit()

            if config.has_option("Flexbb", "path_receptors"):
                self.path_receptors = self.p + config["Flexbb"].get("path_receptors")
            else:
                self.logger.info(
                    "path for receptor flexible BB not found. Use 'path_receptors' parameter"
                )
                exit()

    # --- Input Params -----------------------------+
    def load_refine_docking_parameters(self, config):
        if self.docking_type_option == "Refine":
            if config.has_option("Inputs", "init_pdbs"):
                self.init_pdbs = config["Inputs"].get("init_pdbs")
            else:
                self.logger.info(
                    "Models used for refinement are not found. Use 'init_pdbs' parameter"
                )
                exit()

    def load_output_parameters(self, config):
        # --- OUTPUT PATH -----------------------------------+
        self.out_path = config["Outputs"].get("output_path", fallback="")
        self.out_path = "./" + self.out_path
        if self.out_path[-1] != "/":
            self.out_path = self.out_path + "/"
        os.makedirs(self.out_path, exist_ok=True)

        # ---- OUTPUT PDB -----------------------------------+
        self.out_pdb = config["Outputs"].get("output_pdb", fallback="False")
        self.out_pdb = ast.literal_eval(self.out_pdb)
