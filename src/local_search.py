import logging

from pyrosetta import Pose
from src.individual import Individual
from src.local_search_strategy import LocalSearchStrategy
from src.utils import get_position_info

class LocalSearchPopulation:
    def __init__(self, scfxn, packer_option, config, dockmetric):
        self.config = config
        self.scfxn = scfxn
        self.local_logger = logging.getLogger("evodock.local")
        self.local_logger.setLevel(logging.INFO)
        self.local_search_strategy = LocalSearchStrategy(config, scfxn, scfxn.dock_pose)
        # self.starting_pose = Pose()
        # self.starting_pose.assign(scfxn.dock_pose)
        self.best_pose = Pose()
        self.best_pose.assign(scfxn.dock_pose)
        self.best_score = 100000000
        self.dockmetric = dockmetric

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def process_individual(self, ind, local_search=True):
        data = self.local_search_strategy.apply(ind, local_search)
        pose = data["pose"]
        rmsd = self.dockmetric.ca_rmsd(pose)
        interface = self.dockmetric.interaction_energy(pose)
        irms = self.dockmetric.interface_rmsd(pose)
        if data["after"] < self.best_score:
            self.best_score = data["after"]
            self.best_pose.assign(pose)
        # get position from pose
        positions = get_position_info(pose, self.config.syminfo)
        # replace trial with this new positions
        genotype = self.scfxn.convert_positions_to_genotype(positions)
        result_individual = Individual(
        idx=ind.idx,
        genotype= genotype,
        score = data["after"],
        idx_ligand= data["idx_ligand"],
        idx_receptor=data["idx_receptor"],
        idx_subunit=data["idx_subunit"],
        rmsd=rmsd,
        i_sc=interface,
        irms=irms,
        ligand_name=ind.ligand_name,
        receptor_name=ind.receptor_name,
        subunit_name=ind.subunit_name,
        flipped=ind.flipped,
        fixed=ind.fixed
        )
        self.scfxn.dock_pose.assign(self.best_pose)
        return result_individual, data["before"], data["after"]