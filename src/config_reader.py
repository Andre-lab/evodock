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
from pathlib import Path
import pandas as pd
import math

MAIN_PATH = os.getcwd()

class EvodockConfig:
    def __init__(self, ini_file):
        print(ini_file)
        if ".ini" not in ini_file:
            print("config.ini file not found")
            print("commandline: evodock.py <path_to_config.ini>")
            exit()
        self.logger = logging.getLogger("evodock.config")
        self.logger.setLevel(logging.INFO)
        config = self.read_config(ini_file)
        self.load_docking(config)
        self.load_symmetry(config)
        self.load_inputs(config)
        self.load_native(config)
        self.load_outputs(config)
        self.load_bounds(config)
        self.load_flexbb(config)
        self.load_pymol(config)
        self.load_rosetta_options(config)
        self.load_DE(config)
        self.load_refine_docking_parameters(config)
        self.load_seed(config)
        self.config = config

    def load_docking(self, config):
        """Load the docking_type."""
        self.docking_type_option = config["Docking"].get("type")
        if not self.docking_type_option in ("Global, GlobalFromMultimer, Local, Refine"):
            raise ValueError("The docking_type must be of either 'Global', 'GlobalFromMultimer', 'Local' or 'Refine'")
        if self.docking_type_option == "Refine":
            raise NotImplementedError("Refine is not implemented properly as of yet!")
        self.logger.info(f" Docking mode: {self.docking_type_option}")

    def load_inputs(self, config):
        """Loads all information in the Inputs option."""
        self.flexbb = False
        self.path_to_subunit_folder = None
        # check if we are running a single or multiple backbones
        if config.has_option("Inputs", "single"):
            self.pose_input = config["Inputs"].get("single")
            self.logger.info(f" Single mode detected (single={self.pose_input})")
        elif config.has_option("Inputs", "subunits"):
            # find path to either a single structure a multiple structures
            self.path_to_subunit_folder = config["Inputs"].get("subunits")
            self.subunit_paths = glob.glob(self.path_to_subunit_folder + "*")
            self.subunit_library_size = len(self.subunit_paths)
            assert self.subunit_library_size > 0, f"The size of the subunit library is less than 1. At least 1 or more subunit files must be present in {self.path_to_subunit_folder}"
            self.pose_input = random.choice(self.subunit_paths)
            self.logger.info(f" Ensemble mode detected (subunits={self.path_to_subunit_folder}, size={self.subunit_library_size}, initial backbone={Path(self.pose_input).name})")
            self.flexbb = True
        else:
            if self.symmetric:
                self.logger.info(" For a symmetric system. Use either the 'subunits' or the 'single' parameter")
                exit()
            # find path to either a single structure a multiple structures
            if config.has_option("Inputs", "path_ligands") and config.has_option("Inputs", "path_receptors"):
                self.path_to_ligand_folder = config["Inputs"].get("path_ligands")
                self.path_to_receptor_folder =  config["Inputs"].get("path_receptors")
                self.ligand_paths = glob.glob(self.path_to_ligand_folder + "*")
                self.receptor_paths = glob.glob(self.path_to_receptor_folder + "*")
                self.ligand_library_size = len(self.ligand_paths)
                self.receptor_library_size = len(self.ligand_paths)
                assert self.ligand_library_size > 0, "The size of the ligand library is less than 1. At least 1 or more ligand files must be present"
                assert self.receptor_library_size > 0, "The size of the receptor library is less than 1. At least 1 or more receptor files must be present"
                self.flexbb = True
                self.flexbb = True
                raise NotImplementedError("WE SHOULD PICK AN INITIAL POSE INPUT - SO we have to combine a receptor and ligand right?")
                self.pose_input = random.choice(glob.glob(self.path_to_subunit_folder + "*"))
            else:
                self.logger.info(" For a non-symmetric system. Use both 'path_ligands' and 'path_receptors' or the 'single' parameter")
                exit()
        # get the symmetry file
        if self.symmetric:
            self.syminfo.input_symdef = config.get("Inputs", "symdef_file")

    def load_native(self, config):
        """Loads all information in the Native option."""
        self.native_input = None
        if config["Native"].get("crystallic_input"):
            self.native_input = config["Native"].get("crystallic_input")
            if self.symmetric and config.has_option("Native", "symmetric_input"):
                self.syminfo.native_symmetric_input = config["Native"].get("symmetric_input")
                self.syminfo.native_symdef = config.get("Native", "symdef_file")

    def load_bounds(self, config):
        """Loads all information in the Bounds option."""
        if self.symmetric:
            if config.has_option("Bounds", "init"):
                self.syminfo.init_bounds = [(-float(i), float(i)) for i in config.get("Bounds", "init").split(",")]
                # init_bounds = iter([i for i in config.get("Bounds", "init").split(",")])
                # init_bounds = [list(islice(init_bounds, l)) for l in [len(i) for i in self.syminfo.dofs_str]]
                # self.syminfo.init_bounds = [(-l, l) for l in [float(i[0]) / float(b[0]) for b, i in zip(self.syminfo.bounds, init_bounds)]]
            else:
                self.syminfo.init_bounds = None
            self.syminfo.bounds = [float(i) for i in config.get("Bounds", "bounds").split(",")]
            self.syminfo.x_transfile = config.get("Bounds", "xtrans_file", fallback=None)
            if self.syminfo.x_transfile is not None:
                self.syminfo.x_transfile = pd.read_csv(self.syminfo.x_transfile)
            self.syminfo.normalize_trans = [int(i) for i in config.get("Bounds", "normalize_trans", fallback="2000,1000").split(",")]
            self.syminfo.bound_penalty = config.getfloat("Bounds", "bound_penalty", fallback=1)
            self.allow_flip = config.getboolean("Bounds", "allow_flip", fallback=False)
            self.init_input_fix_percent = config.getfloat("Bounds", "init_input_fix_percent", fallback=0) / 100
            assert self.init_input_fix_percent >= 0 or self.init_input_fix_percent <= 1, "init_input_fix_percent has to be between 0 and 100."
            if math.isclose(self.init_input_fix_percent, 0):
                self.init_input_fix_percent = None
        else:
            raise NotImplementedError

    def load_flexbb(self, config):
        if config.has_section("Flexbb"):
            self.swap_prob = config["Flexbb"].getfloat("swap_prob", fallback=0.3)
            self.low_memory_mode = config["Flexbb"].getboolean("low_memory_mode", fallback=True)
        else:
            self.swap_prob = None
            self.low_memory_mode = False


    def load_outputs(self, config):
        self.out_path = config["Outputs"].get("output_path", fallback="")
        os.makedirs(self.out_path, exist_ok=True)
        self.out_pdb = config["Outputs"].get("output_pdb", fallback=True)
        self.out_pdb = ast.literal_eval(self.out_pdb)
        self.output_pdb_per_generation = config["Outputs"].getboolean("output_pdb_per_generation", fallback=False)

    def load_rosetta_options(self, config):
        """Loads in Rosetta options"""
        self.rosetta_options = None
        if config.has_section("RosettaOptions"):
            # this is fed into init later on
            self.rosetta_options = config.items("RosettaOptions")

    def load_seed(self, config):
        """Loads the seed if set in the config file."""
        if config.has_option("Seed", "seed"):
            self.seed = config.getint("Seed", "seed")
        else:
            self.seed = None

    # -------- visualizers  -------- #

    def visualize_pose(self, pose, idx, extra=None):
        """Visualizes the pose in PyMOL."""
        if self.pmm:
            name = str(idx)
            name = name + f"_{extra}" if extra else name
            pose.pdb_info().name(name)
            self.pmm.apply(pose)

    def visualize_genotype(self, genotype, scfxn):
        """Visualize genotype. Only here for testing purposes."""
        if self.pmm is None:
            self.pmm = PyMOLMover(address=self.ipaddress, port=65000, max_packet_size=1400)
        scfxn.native_pose.pdb_info().name("native")
        self.pmm.apply(scfxn.native_pose)
        pose = scfxn.apply_genotype_to_pose(genotype)
        name = "pose_from_genotype"
        self.pmm.apply(pose)
        ...

    def visualize_population(self, population, scfxn):
        """Visualize population. Only here for testing purposes."""
        if self.pmm is None:
            self.pmm = PyMOLMover(address=self.ipaddress, port=65000, max_packet_size=1400)
        scfxn.native_pose.pdb_info().name("native")
        self.pmm.apply(scfxn.native_pose)
        for ind in population:
            pose = scfxn.apply_genotype_to_pose(ind.genotype)
            name = str(ind.idx)
            pose.pdb_info().name(name)
            self.pmm.apply(pose)

    def load_symmetry(self, config):
        """Checks for symmetry and constructs a SymInfo object containing all symmetry information."""
        if config.has_option("Inputs", "symdef_file"):
            self.logger.info(" Symmetry is detected.")
            self.symmetric = True
            self.syminfo = SymInfo()
        else:
            self.syminfo = None
            self.symmetric = False
            self.logger.info(" Symmetry not detected.")

    def load_pymol(self, config):
        self.ipaddress = config.get("Pymol", "ipaddress", fallback="0.0.0.0")
        if config.has_option("Pymol", "on") and config.getboolean("Pymol", "on"):
            self.pmm = PyMOLMover(address=self.ipaddress, port=65000, max_packet_size=1400)
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

    def load_DE(self, config):
        # --- DE PARAMS -----------------------------------+
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
        # --- SLIDE PARAMS -----------------------------------+
        if config.has_option("DE", "slide"):
            self.slide = config.getboolean("DE", "slide")
            self.max_slide_attempts = config["DE"].getfloat("max_slide_attempts", fallback=100)
            self.slide_trans_mag = config["DE"].getfloat("slide_trans_mag", fallback=0.3)
        else:
            self.slide = True

    def load_refine_docking_parameters(self, config):
        if self.docking_type_option == "Refine":
            if config.has_option("Refine", "init_pdbs"):
                self.init_pdbs = self.p + config["Refine"].get("init_pdbs")
            else:
                self.logger.info(
                    "Models used for refinement are not found. Use 'init_pdbs' parameter"
                )
                exit()


