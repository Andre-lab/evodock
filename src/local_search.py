import glob
import logging

from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.core.import_pose import poses_from_files
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.protocols.docking import (DockMCMProtocol,
                                                 FaDockingSlideIntoContact,
                                                 calc_interaction_energy,
                                                 calc_Irmsd)
from pyrosetta.rosetta.utility import vector1_std_string

from src.individual import Individual
from src.utils import get_position_info


class LocalSearchPopulation:
    # Options:
    # None: only score and return the poses
    # only_slide: just slide_into_contact
    # mcm_rosetta: mcm protocol mover (high res) from rosetta (2 cycles)
    def __init__(self, scfxn, packer_option="default_combination"):
        lst_ligand = glob.glob(
            "/home/daniel/projects/evodock/flexbackbones/relax/2jtoA/*"
        )
        lst_ligand += lst_ligand
        lst_ligand += lst_ligand
        lst_ligand += lst_ligand
        lst_ligand = lst_ligand[:3]

        lst_receptor = glob.glob(
            "/home/daniel/projects/evodock/flexbackbones/relax/1kwmA/*"
        )
        lst_receptor += lst_receptor
        lst_receptor += lst_receptor
        lst_receptor += lst_receptor
        lst_receptor = lst_receptor[:3]

        filenames_ligand = vector1_std_string()
        for f in lst_ligand:
            filenames_ligand.append(f)
        filenames_receptor = vector1_std_string()
        for f in lst_receptor:
            filenames_receptor.append(f)

        self.list_ligand = poses_from_files(filenames_ligand)
        self.list_receptor = poses_from_files(filenames_receptor)
        self.packer_option = packer_option
        self.scfxn = scfxn
        self.local_logger = logging.getLogger("evodock.local")
        self.local_logger.setLevel(logging.INFO)
        self.native_fold_tree = scfxn.dock_pose.fold_tree()
        # self.slide_into_contact = DockingSlideIntoContact(1)
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
        # pose2 = self.list_ligand[1]
        # pose1 = self.list_receptor[1]
        # join_pose = Pose()
        # join_pose.assign(pose1)
        # append_pose_to_pose(join_pose, pose2, True)
        # # print("before : ")
        # # print(join_pose.fold_tree())
        # join_pose.conformation().detect_disulfides()
        # # print("after : ")
        # # print(join_pose.fold_tree())

        # join_pose.fold_tree(self.native_fold_tree)
        # self.scfxn.dock_pose = join_pose

        pose = self.scfxn.apply_genotype_to_pose(ind)

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
            pose, self.scfxn.scfxn_rosetta, Vector1([1])
        )
        irms = calc_Irmsd(
            self.scfxn.native_pose, pose, self.scfxn.scfxn_rosetta, Vector1([1])
        )

        # get position from pose
        positions = get_position_info(pose)
        # replace trial with this new positions
        genotype = self.scfxn.convert_positions_to_genotype(positions)
        result_individual = Individual(genotype, after, 1, 1, rmsd, interface, irms)
        return result_individual, before, after
