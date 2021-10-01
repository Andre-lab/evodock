import glob
import logging
import random
from random import sample

from pyrosetta import Pose, Vector1
from pyrosetta.rosetta.core.import_pose import poses_from_files
from pyrosetta.rosetta.core.pose import (addVirtualResAsRoot,
                                         append_pose_to_pose, chain_end_res,
                                         remove_virtual_residues)
from pyrosetta.rosetta.core.scoring import CA_rmsd, calpha_superimpose_pose
from pyrosetta.rosetta.protocols.docking import (DockMCMProtocol,
                                                 FaDockingSlideIntoContact,
                                                 calc_interaction_energy,
                                                 calc_Irmsd, setup_foldtree)
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.utility import vector1_std_string

from src.individual import Individual
from src.utils import get_pose_from_file, get_position_info


def apply_bb_strategy(bb_strategy, local_search, ind, pose):
    if bb_strategy == "library":
        join_pose, idx_receptor, idx_ligand = local_search.define_ensemble(ind, pose)
    if bb_strategy == "relax":
        join_pose, idx_receptor, idx_ligand = local_search.define_relaxedbackbone(pose)
    if bb_strategy == "fixed":
        join_pose = pose
        idx_receptor, idx_ligand = 1, 1
    return join_pose, idx_receptor, idx_ligand


class LocalSearchPopulation:
    # Options:
    # None: only score and return the poses
    # only_slide: just slide_into_contact
    # mcm_rosetta: mcm protocol mover (high res) from rosetta (2 cycles)
    def __init__(self, scfxn, packer_option, config):
        self.config = config
        lst_ligand = glob.glob(config.path_ligands)
        # lst_ligand = lst_ligand[:10]

        lst_receptor = glob.glob(config.path_receptors)
        # lst_receptor = lst_receptor[:10]

        filenames_ligand = vector1_std_string()
        for f in lst_ligand:
            filenames_ligand.append(f)
        filenames_receptor = vector1_std_string()
        for f in lst_receptor:
            filenames_receptor.append(f)

        self.list_ligand = poses_from_files(filenames_ligand)
        self.list_receptor = poses_from_files(filenames_receptor)
        self.ref_ligand = get_pose_from_file("./easy_dock/ligand.kb.clean.pdb")
        self.ref_receptor = get_pose_from_file("./easy_dock/receptor.clean.kb.pdb")
        print("list ligand {}".format(len(self.list_ligand)))
        print("list receptor {}".format(len(self.list_receptor)))
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

        # self.docksetup = PrepareForDocking(scfxn.dock_pose)
        # pose = PrepareForDocking(scfxn.dock_pose).apply(scfxn.dock_pose)
        self.native_fold_tree = scfxn.dock_pose.fold_tree()
        self.starting_pose = Pose()
        self.starting_pose.assign(scfxn.dock_pose)
        self.best_pose = Pose()
        self.best_pose.assign(scfxn.dock_pose)
        self.best_score = 100000000
        self.relax = FastRelax(0)
        self.relax.set_scorefxn(self.scfxn.scfxn_rosetta)
        self.relax.max_iter(1)
        self.relaxed_backbones = [self.scfxn.dock_pose]
        self.bb_strategy = config.bb_strategy

    def energy_score(self, pose):
        score = self.scfxn.scfxn_rosetta(pose)
        return score

    def define_ensemble(self, ind, reference_pose):
        if random.uniform(0, 1) < 0.5:
            idx_receptor = ind.idx_receptor
            idx_ligand = random.randint(1, len(self.list_ligand))
        else:
            idx_receptor = random.randint(1, len(self.list_receptor))
            idx_ligand = ind.idx_ligand
        pose1 = self.list_receptor[idx_receptor]
        pose2 = self.list_ligand[idx_ligand]

        pose_R = Pose(reference_pose, 1, chain_end_res(reference_pose, 1))
        pose_L = Pose(
            reference_pose,
            chain_end_res(reference_pose, 1) + 1,
            reference_pose.total_residue(),
        )
        pose_R.pdb_info().name("chainA")
        pose_L.pdb_info().name("chainB")

        calpha_superimpose_pose(pose1, pose_R)
        calpha_superimpose_pose(pose2, pose_L)

        join_pose = Pose()
        join_pose.assign(pose1)
        append_pose_to_pose(join_pose, pose2, True)

        # addVirtualResAsRoot(join_pose)

        # self.docksetup.apply(join_pose)
        join_pose.fold_tree(self.native_fold_tree)
        join_pose.conformation().detect_disulfides()
        return join_pose, idx_receptor, idx_ligand

    def define_relaxedbackbone(self, pose):
        idx_receptor, idx_ligand = 1, 1
        join_pose = self.relaxed_backbones[
            random.randint(0, len(self.relaxed_backbones) - 1)
        ]

        calpha_superimpose_pose(join_pose, pose)
        join_pose.fold_tree(self.native_fold_tree)
        join_pose.conformation().detect_disulfides()
        return join_pose, idx_receptor, idx_ligand

    def process_individual(self, ind, local_search=True):
        pose = self.scfxn.apply_genotype_to_pose(ind.genotype)
        if self.bb_strategy == "library":
            join_pose, idx_receptor, idx_ligand = apply_bb_strategy(
                self.bb_strategy, self, ind, pose
            )
        else:
            rand = random.uniform(0, 1)
            if rand > 0.1:
                join_pose, idx_receptor, idx_ligand = apply_bb_strategy(
                    self.bb_strategy, self, ind, pose
                )
            else:
                join_pose, idx_receptor, idx_ligand = apply_bb_strategy(
                    self.bb_strategy, self, ind, pose
                )

        self.scfxn.dock_pose = join_pose
        slide_into_contact = FaDockingSlideIntoContact(pose.num_jump())
        # pose.dump_pdb("test_join.pdb")
        before = self.energy_score(pose)
        if local_search and self.packer_option != "None":
            # slide_into_contact.apply(pose)
            if self.bb_strategy == "library":
                rand = 1.0
            else:
                rand = random.uniform(0, 1)
            if rand > 0.1:
                self.docking.apply(pose)
            else:
                self.relax.apply(pose)
                self.relaxed_backbones.append(pose)
                if len(self.relaxed_backbones) > 100:
                    self.relaxed_backbones = sample(self.relaxed_backbones, 100)
                remove_virtual_residues(pose)
                pose.fold_tree(self.native_fold_tree)
            after = self.energy_score(pose)

        else:
            after = before

        rmsd = self.scfxn.get_rmsd(pose)

        interface = calc_interaction_energy(
            pose, self.scfxn.scfxn_rosetta, Vector1([1])
        )
        irms = calc_Irmsd(
            self.scfxn.native_pose, pose, self.scfxn.scfxn_rosetta, Vector1([1]),
        )

        if after < self.best_score:
            self.best_score = after
            self.best_pose.assign(pose)

        # get position from pose
        positions = get_position_info(pose)
        # replace trial with this new positions
        genotype = self.scfxn.convert_positions_to_genotype(positions)

        result_individual = Individual(
            genotype, after, idx_ligand, idx_receptor, rmsd, interface, irms
        )

        self.scfxn.dock_pose.assign(self.starting_pose)
        return result_individual, before, after
