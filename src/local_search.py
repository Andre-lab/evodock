import glob
import logging
import random

from pyrosetta import Pose, SwitchResidueTypeSetMover, Vector1
from pyrosetta.rosetta.core.import_pose import poses_from_files
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.docking import (ConformerSwitchMover,
                                                 DockingEnsemble,
                                                 DockMCMProtocol,
                                                 FaDockingSlideIntoContact,
                                                 calc_interaction_energy,
                                                 calc_Irmsd)
from pyrosetta.rosetta.protocols.simple_moves import ReturnSidechainMover
from pyrosetta.rosetta.utility import vector1_std_string

from src.individual import Individual
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
    # Options:
    # None: only score and return the poses
    # only_slide: just slide_into_contact
    # mcm_rosetta: mcm protocol mover (high res) from rosetta (2 cycles)
    def __init__(self, scfxn, packer_option, config):
        self.config = config
        self.conformer1_, self.conformer2_ = ConformerSwitchCreator(config)
        self.packer_option = packer_option
        self.scfxn = scfxn
        self.local_logger = logging.getLogger("evodock.local")
        self.local_logger.setLevel(logging.INFO)
        self.native_fold_tree = scfxn.dock_pose.fold_tree()
        # self.slide_into_contact = DockingSlideIntoContact(1)

        self.recover_sidechains = ReturnSidechainMover(scfxn.native_pose)
        self.to_centroid = SwitchResidueTypeSetMover("centroid")
        self.to_fullatom = SwitchResidueTypeSetMover("fa_standard")
        if packer_option == "mcm_rosetta":
            mcm_docking = DockMCMProtocol()
            mcm_docking.set_native_pose(scfxn.dock_pose)
            mcm_docking.set_scorefxn(scfxn.scfxn_rosetta)
            mcm_docking.set_rt_min(False)
            mcm_docking.set_sc_min(False)
            mock_pose = Pose()
            mock_pose.assign(scfxn.dock_pose)
            mcm_docking.apply(mock_pose)
            self.docking = mcm_docking
            self.docking.set_task_factory(mcm_docking.task_factory())
            self.docking.set_ignore_default_task(True)

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
        idx_ligand, idx_receptor, join_pose = self.build_essemble_pose(ind)
        self.scfxn.dock_pose = join_pose
        # idx_ligand, idx_receptor = 1, 1

        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)

        slide_into_contact = FaDockingSlideIntoContact(pose.num_jump())
        pose.dump_pdb("test_join.pdb")
        before = self.energy_score(pose)
        if local_search and self.packer_option != "None":
            slide_into_contact.apply(pose)
            if self.packer_option != "only_slide":
                self.docking.apply(pose)
            after = self.energy_score(pose)
        else:
            after = before

        rmsd = self.scfxn.get_rmsd(pose)

        interface = calc_interaction_energy(
            pose, self.scfxn.scfxn_rosetta, Vector1([pose.num_jump()])
        )
        irms = calc_Irmsd(
            self.scfxn.native_pose,
            pose,
            self.scfxn.scfxn_rosetta,
            Vector1([pose.num_jump()]),
        )

        if after < self.best_score:
            self.best_pose.assign(pose)

        # get position from pose
        positions = get_position_info(pose)
        # replace trial with this new positions
        genotype = self.scfxn.convert_positions_to_genotype(positions)
        result_individual = Individual(
            genotype, after, idx_ligand, idx_receptor, rmsd, interface, irms
        )

        return result_individual, before, after
