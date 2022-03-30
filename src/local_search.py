import logging

from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.protocols.docking import calc_interaction_energy, calc_Irmsd

from src.individual import Individual
from src.local_search_strategy import LocalSearchStrategy
from src.utils import get_position_info


class LocalSearchPopulation:
    def __init__(self, scfxn, packer_option, config):
        self.config = config
        self.scfxn = scfxn
        self.local_logger = logging.getLogger("evodock.local")
        self.local_logger.setLevel(logging.INFO)
        self.local_search_strategy = LocalSearchStrategy(config, scfxn, scfxn.dock_pose)
        self.starting_pose = Pose()
        self.starting_pose.assign(scfxn.dock_pose)
        self.best_pose = Pose()
        self.best_pose.assign(scfxn.dock_pose)
        self.best_score = 100000000
        self.best_pose = Pose()
        self.best_pose.assign(scfxn.dock_pose)
        self.best_score = 1000

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def process_individual(self, ind, local_search=True):
        data = self.local_search_strategy.apply(ind, local_search)
        pose = data["pose"]
        rmsd = self.scfxn.get_rmsd(pose)
        interface = calc_interaction_energy(
            pose, self.scfxn.scfxn_rosetta, Vector1([1])
        )
        irms = calc_Irmsd(
            self.scfxn.native_pose,
            pose,
            self.scfxn.scfxn_rosetta,
            Vector1([1]),
        )
        if data["after"] < self.best_score:
            self.best_score = data["after"]
            self.best_pose.assign(pose)
        # get position from pose
        positions = get_position_info(pose)
        # replace trial with this new positions
        genotype = self.scfxn.convert_positions_to_genotype(positions)
        result_individual = Individual(
            genotype,
            data["after"],
            data["idx_ligand"],
            data["idx_receptor"],
            rmsd,
            interface,
            irms,
        )
        # if self.config.bb_strategy == "only_relax":
        #     if self.scfxn.score(self.scfxn.dock_pose) > data["after"]:
        #         self.scfxn.dock_pose.assign(pose)
        # else:

        self.scfxn.dock_pose.assign(self.starting_pose)

        return result_individual, data["before"], data["after"]
