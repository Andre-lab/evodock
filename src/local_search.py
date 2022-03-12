import logging

from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.protocols.docking import calc_interaction_energy, calc_Irmsd

from src.individual import Individual
from src.local_search_strategy import LocalSearchStrategy
from src.utils import get_position_info


# start_res_ 1 end_res_ 306 conf_size_ 306
# start_res_ 307 end_res_ 381 conf_size_ 75 ensemble_size_ 101
def ConformerSwitchCreator(config):
    scorelow = ScoreFunctionFactory.create_score_function("score3")
    scorehigh = ScoreFunctionFactory.create_score_function("ref2015")
    ensemble1_ = DockingEnsemble(
        1, 306, 1, config.path_receptors, "A", scorelow, scorehigh
    )
    ensemble2_ = DockingEnsemble(
        307, 381, 1, config.path_ligands, "B", scorelow, scorehigh
    )

    conformer1_ = ConformerSwitchMover(ensemble1_, True)
    conformer2_ = ConformerSwitchMover(ensemble2_, True)
    return conformer1_, conformer2_


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

    def build_essemble_pose(self, ind):
        if random.uniform(0, 1):
            idx_receptor = ind.idx_receptor
            idx_ligand = random.randint(1, 100)
        else:
            idx_ligand = ind.idx_ligand
            idx_receptor = random.randint(1, 100)

        join_pose = Pose()
        join_pose.assign(self.scfxn.input_pose)
        self.to_centroid.apply(join_pose)
        self.conformer1_.switch_conformer(join_pose, idx_receptor)
        self.conformer2_.switch_conformer(join_pose, idx_ligand)
        self.to_fullatom.apply(join_pose)
        self.recover_sidechains.apply(join_pose)
        return idx_ligand, idx_receptor, join_pose

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
        #     if score(self.scfxn.dock_pose) > data["after"]:
        #         self.scfxn.dock_pose.assign(pose)
        # else:

        self.scfxn.dock_pose.assign(self.starting_pose)

        return result_individual, data["before"], data["after"]
