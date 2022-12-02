#!/usr/bin/env python
# coding: utf-8
import configparser
import logging
import os
import ast
from src.utils import get_translation_max
from itertools import islice
from pyrosetta import PyMOLMover
from src.individual import Individual
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
import glob
import random
from src.symmetry import SymInfo


MAIN_PATH = os.getcwd()

class EvodockConfig:
    def __init__(self, ini_file):
        print(ini_file)
        if ".ini" not in ini_file:
            print("config.ini file not found")
            print("commandline: evodock.py <path_to_config.ini>")
            exit()

        self.p = os.getcwd()
        self.logger = logging.getLogger("evodock.config")
        self.logger.setLevel(logging.INFO)


        config = self.read_config(ini_file)

        # depends on symmetry
        self.load_symmetry_parameters(config)

        # depend on docking type
        self.load_flexible_docking_parameters(config)
        self.load_refine_docking_parameters(config)

        # FIXME THIS NEEDS TO BE REWRITTEN
        self.check_required_parameters(config)
        self.load_required_parameters(config)
        self.load_output_parameters(config)
        self.load_seed(config)


        # general docking parameters
        self.load_docking_parameters(config)

        # minimization options
        self.load_minimization_options(config)

        # depends on pymol
        self.load_pymol_parameters(config)

        # --- CONFIG STRUCTURE -----------------------------------+
        self.config = config

    def load_minimization_options(self, config):
        self.minimization_option = config.get("Minimization", "minimization_option", fallback=None)
        self.min_option_bb =  config.get("Minimization", "min_option_bb", fallback=None)
        self.min_option_sc = config.getboolean("Minimization", "min_option_sc", fallback=False)
        self.cartesian = config.getboolean("Minimization", "cartesian", fallback=False)

    def load_seed(self, config):
        """Loads the seed if set in the config file."""
        if config.has_option("Seed", "seed"):
            self.seed = config.getint("Seed", "seed")
        else:
            self.seed = None

    def visualize_pose(self, pose, idx, extra=None):
        """Visualizes the pose in PyMOL."""
        if self.pmm:
            name = str(idx)
            name = name + f"_{extra}" if extra else name
            pose.pdb_info().name(name)
            self.pmm.apply(pose)

    def load_docking_parameters(self, config):
        """Load docking parameters"""
        if config.has_option("Symmetry", "input_symdef_file"):
            self.num_first_cycle = config.getint("DE", "num_first_cycle", fallback=1)
            self.num_second_cycle = config.getint("DE", "num_second_cycle", fallback=1)
        else:
            self.num_first_cycle = config.getint("DE", "num_first_cycle", fallback=1)
            self.num_second_cycle = config.getint("DE", "num_second_cycle", fallback=1)

    def load_symmetry_parameters(self, config):
        """Extracts the symmetrical information from the config file."""
        if config.has_option("Symmetry", "input_symdef_file"):
            self.syminfo = SymInfo()
            self.syminfo.store_info_from_config(config)
        else:
            self.syminfo = None

    def load_pymol_parameters(self, config):
        if config.has_option("Pymol", "on") and config.getboolean("Pymol", "on"):
            self.ipaddress = config.get("Pymol", "ipaddress", fallback="0.0.0.0")
            self.pmm = PyMOLMover(address=self.ipaddress, port=65000, max_packet_size=1400)
            if config.has_option("Pymol", "history"):
                self.pmm.keep_history(config.getboolean("Pymol", "history"))
        else:
            self.pmm = None
            self.ipaddress = None
        if config.has_option("Pymol", "show_local_search"):
            self.show_local_search = config.getboolean("Pymol", "show_local_search")
        else:
            self.show_local_search = False

    def get_max_translation(self):
        # --- Position Params -----------------------------+
        if not self.syminfo:
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
            ("Inputs", "native_input"),
            ("DE", "scheme"),
            ("DE", "popsize"),
            ("DE", "mutate"),
            ("DE", "recombination"),
            ("DE", "maxiter"),
            ("DE", "local_search"),
        ]
        if self.docking_type_option != "Flexbb":
            req_parameters.append(("Inputs", "pose_input"))
        all_required_found = True
        # for param in req_parameters:
        #     if config.has_option(param[0], param[1]):
        #         self.docking_type_option = config[param[0]].get(param[1])
        #     else:
        #         self.logger.info(f"[{param[0]}]{param[1]} not found.")
        #         all_required_found = False
        if not all_required_found:
            exit()

    def load_required_parameters(self, config):

        # -- INPUT PARAMETERS ---------------------------------------+
        # load a random subunit/ligand/receptor if using flexbb
        if self.docking_type_option != "Flexbb":
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

        # symmetry specific
        if self.syminfo:
            native_input_symmetric_file_abs = config["Inputs"].get("native_symmetric_input")
            native_input_symmetric_file_relative = self.p + "/" + config["Inputs"].get("native_symmetric_input")
            if os.path.isfile(native_input_symmetric_file_abs):
                self.native_symmetric_input = native_input_symmetric_file_abs
            else:
                self.native_symmetric_input = native_input_symmetric_file_relative

            if not os.path.isfile(self.native_symmetric_input):
                self.logger.info(f"native file not found: {self.native_symmetric_input}")
                exit()

        # --- DE PARAMS ---------------------------------------------+
        self.scheme = config["DE"].get("scheme")
        self.popsize = config["DE"].getint("popsize")
        self.mutate = config["DE"].getfloat("mutate")
        self.recombination = config["DE"].getfloat("recombination")
        self.maxiter = config["DE"].getint("maxiter")

        if config.has_option("DE", "selection"):
            self.selection = config["DE"].get("selection")
            assert self.selection in ("total", "interface"), "only 'total' and 'interface' is understood for DE:selection."
        else:
            self.selection = "total"


        # --- MEMETIC PARAMS -----------------------------------+
        if config.has_option("DE", "local_search"):
            self.local_search_option = config["DE"].get("local_search")
        else:
            self.logger.info("DANGER: local_search is None")
            self.local_search_option = "None"

        # --- SLIDE option -----------------------------------+
        if config.has_option("DE", "slide"):
            self.slide = config.getboolean("DE", "slide")
            self.max_slide_attempts = config["DE"].getboolean("max_slide_attempts", fallback=100)
            self.slide_trans_mag = config["DE"].getboolean("slide_trans_mag", fallback=0.3)
        else:
            self.slide = True

    def load_flexible_docking_parameters(self, config):
        self.docking_type_option = config["Docking"].get("type")
        if self.docking_type_option == "Flexbb":
            self.swap_prob = config["Flexbb"].getfloat("swap_prob", fallback=0.3)
            self.low_memory_mode = config["Flexbb"].getboolean("low_memory_mode", fallback=False)
            self.normalize_score = config["Flexbb"].getboolean("normalize_score", fallback=False)
            self.save_bbs = config["Flexbb"].get("save_bbs", fallback=None)
            self.dssp_acceptance = config["Flexbb"].getfloat("dssp_acceptance", fallback=95.0)
            self.rmsd_acceptance = config["Flexbb"].getfloat("rmsd_acceptance", fallback=2.0)
            if self.syminfo:
                self.path_subunits = self.p + config["Flexbb"].get("subunits")
                self.pose_input = random.choice(glob.glob(self.path_subunits + "*"))
                self.logger.info(
                    f"Selected {self.pose_input} as the starting structure."
                )
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
        else:
            self.swap_prob = None
            self.low_memory_mode = None
            self.normalize_score = None
            self.save_bbs = None
            self.dssp_acceptance = None
            self.rmsd_acceptance = None
            self.normalize_score = False

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
        self.out_pdb = config["Outputs"].get("output_pdb", fallback=True)
        self.out_pdb = ast.literal_eval(self.out_pdb)
        self.output_pdb_per_generation = config["Outputs"].getboolean("output_pdb_per_generation", fallback=False)
