#!/usr/bin/env python
# coding: utf-8
import configparser
import logging
import os
import ast
from src.utils import get_translation_max
from itertools import islice
from pyrosetta import PyMOLMover
from src.utils import IP_ADDRESS

MAIN_PATH = os.getcwd()

class EvodockConfig:
    def __init__(self, ini_file):
        if ".ini" not in ini_file:
            print("config .ini file not found")
            print("commandline: evodock.py <path_to_config.ini>")
            exit()

        self.p = os.getcwd()
        self.logger = logging.getLogger("evodock.config")
        self.logger.setLevel(logging.INFO)

        config = self.read_config(ini_file)

        self.check_required_parameters(config)
        self.load_required_parameters(config)
        self.load_output_parameters(config)

        # depends on symmetry
        self.load_symmetry_parameters(config)

        # depend on docking type
        self.load_flexible_docking_parameters(config)
        self.load_refine_docking_parameters(config)

        # depends on pymol
        self.load_pymol_parameters(config)

        # --- CONFIG STRUCTURE -----------------------------------+
        self.config = config

    def load_symmetry_parameters(self, config):
        """Extracts the symmetrical information from the config file."""
        syminfo = {}
        if config.has_option("Symmetry", "input_symdef_file"):
            syminfo["input_symdef"] = MAIN_PATH + config.get("Symmetry", "input_symdef_file")
            syminfo["native_symdef"] = MAIN_PATH + config.get("Symmetry", "native_symdef_file")
            if config.has_option("Symmetry", "symdofs"):
                syminfo["jumps_str"] = [i.split(":")[0] for i in config.get("Symmetry", "symdofs").split(",")]
                syminfo["dofs_str"] = [i.split(":")[1:] for i in config.get("Symmetry", "symdofs").split(",")]
                bounds = iter([i for i in config.get("Symmetry", "symbounds").split(",")])
                # for nicer mapping we have to map the bounds the same ways as the dofs
                syminfo["bounds"] = [list(islice(bounds, l)) for l in [len(i) for i in syminfo["dofs_str"]]]
                syminfo["genotype_size"] = config.get("Symmetry", "symdofs").count(":")
                assert len(syminfo["jumps_str"]) == len(syminfo["dofs_str"])
                assert len(syminfo["jumps_str"]) == len(syminfo["bounds"])
            else: # apply defaults
                raise NotImplementedError
            if config.has_option("Symmetry", "initialize_rigid_body_dofs"):
                syminfo["initialize_rigid_body_dofs"] = config.getboolean("Symmetry", "initialize_rigid_body_dofs")
            else:
                raise NotImplementedError
        self.syminfo = syminfo if syminfo else None # Return None if empty

    def load_pymol_parameters(self, config):
        if config.has_option("Pymol", "on") and config.getboolean("Pymol", "on"):
            self.pmm = PyMOLMover(address=IP_ADDRESS, port=65000, max_packet_size=1400)
            if config.has_option("Pymol", "history"):
                self.pmm.keep_history(config.getboolean("Pymol", "history"))
        else:
            self.pmm = None
        if config.has_option("Pymol", "show_local_search"):
            self.show_local_search = config.getboolean("Pymol", "show_local_search")
        else:
            self.show_local_search = False

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
        pose_input_file_abs = config["Inputs"].get("pose_input")
        pose_input_file_relative = self.p + "/" + config["Inputs"].get("pose_input")
        if os.path.isfile(pose_input_file_abs):
            self.pose_input = pose_input_file_abs
        else:
            self.pose_input = pose_input_file_relative

        if not os.path.isfile(self.pose_input):
            self.logger.info(f"input file not found: {self.pose_input}")
            exit()

        native_input_file_abs = config["Inputs"].get("native_input")
        native_input_file_relative = self.p + "/" + config["Inputs"].get("native_input")
        if os.path.isfile(native_input_file_abs):
            self.native_input = native_input_file_abs
        else:
            self.native_input = native_input_file_relative

        if not os.path.isfile(self.native_input):
            self.logger.info(f"native file not found: {self.native_input}")
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

        # --- SLIDE option -----------------------------------+
        if config.has_option("DE", "slide"):
            self.slide = config.getboolean("DE", "slide")
        else:
            self.slide = True

    def load_flexible_docking_parameters(self, config):
        if self.docking_type_option == "Flexbb":
            if self.syminfo:
                self.path_subunits = self.p + config["Flexbb"].get("subunits")
            else:
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
            if config.has_option("Refine", "init_pdbs"):
                self.init_pdbs = self.p + config["Refine"].get("init_pdbs")
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
        self.output_pdb_per_generation = config["Outputs"].getboolean("output_pdb", fallback="False")
