import logging

from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.protocols.docking import (DockingSlideIntoContact,
                                                 DockMCMProtocol,
                                                 calc_interaction_energy,
                                                 calc_Irmsd)

from src.individual import Individual
from src.utils import get_position_info


class LocalSearchPopulation:
    # Options:
    # None: only score and return the poses
    # only_slide: just slide_into_contact
    # mcm_rosetta: mcm protocol mover (high res) from rosetta (2 cycles)
    def __init__(self, scfxn, packer_option="default_combination"):
        self.packer_option = packer_option
        self.scfxn = scfxn
        self.local_logger = logging.getLogger("evodock.local")
        self.local_logger.setLevel(logging.INFO)
        self.slide_into_contact = DockingSlideIntoContact(1)
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

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def process_individual(self, ind, local_search=True):
        pose = self.scfxn.apply_genotype_to_pose(ind)
        before = self.energy_score(pose)
        if local_search and self.packer_option != "None":
            self.slide_into_contact.apply(pose)
            if self.packer_option != "only_slide":
                self.docking.apply(pose)
            after = self.energy_score(pose)
        else:
            after = before

        rmsd = self.scfxn.get_rmsd(pose)

        interface = calc_interaction_energy(
            pose, self.scfxn.scfxn_rosetta, Vector1([1])
        )
        irms = calc_Irmsd(
            self.scfxn.native_pose, pose, self.scfxn.scfxn_rosetta, Vector1([1])
        )

        # get position from pose
        positions = get_position_info(pose)
        # replace trial with this new positions
        genotype = self.scfxn.convert_positions_to_genotype(positions)
        result_individual = Individual(genotype, after, rmsd, interface, irms)
        return result_individual, before, after
