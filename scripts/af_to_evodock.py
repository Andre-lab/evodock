#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Aligns the alphafold predictions
@Author: Mads Jeppesen
@Date: 5/25/22
"""
import random
from pyrosetta.rosetta.core.pose import residue_center_of_mass, center_of_mass
import numpy as np
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pathlib import Path
from pyrosetta import init, pose_from_file
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.core.scoring import CA_rmsd, CA_rmsd_symmetric
import pandas as pd
import json
import pickle
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from scipy.spatial.distance import cdist
import fnmatch
from cubicsym.alphafold.symmetrymapper import SymmetryMapper
from pyrosetta.rosetta.protocols.loops import get_fa_scorefxn
from prepacking import PosePacker, opts
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from cubicsym.paths import DATA
from cubicsym.cubicsetup import CubicSetup
from cubicsym.actors.cubicsymmetryslider import InitCubicSymmetrySlider
from cloudcontactscore.cloudcontactscorecontainer import CloudContactScoreContainer
from symmetryhandler.reference_kinematics import set_jumpdof_str_str
from cubicsym.alignment import tmalign, contactmap_with_mapped_resi
from distutils.util import strtobool

class AF2EvoDOCK:

    def __init__(self, symmetry, af_output_path, out_dir, plddt_confidence, ss_confidence, disconnect_confidence, avg_plddt_threshold, iptm_plus_ptm_threshold, rmsd_threshold,
                 ensemble, max_monomers=None, max_multimers=None, total_max_models=None, modify_rmsd_to_reach_min_models=None, direction=False, model_type="ranked",
                 use_only_symmetric=True):
        self.symmetry = symmetry
        self.af_output_path = af_output_path
        self.out_dir = out_dir
        self.plddt_confidence = plddt_confidence
        self.ss_confidence = ss_confidence
        self.disconnect_confidence = disconnect_confidence
        self.avg_plddt_threshold = avg_plddt_threshold
        self.iptm_plus_ptm_threshold = iptm_plus_ptm_threshold
        self.init_rmsd_threshold = rmsd_threshold
        self.ensemble = ensemble
        self.max_monomers = max_monomers
        self.max_multimers = max_multimers
        self.total_max_models = total_max_models
        self.prepack = True
        self.modify_rmsd_to_reach_min_models = modify_rmsd_to_reach_min_models
        self.direction = direction
        self.model_type = model_type
        self.use_only_symmetric = use_only_symmetric

    def cut_all_but_chain(self, pose, chain):
        """Cut all chains from the pose except the chain parsed."""
        resi_in_chain = []
        # collect all residues in the chain of the pose
        for resi in range(1, pose.size() + 1):
            if pose.pdb_info().chain(resi) == chain:
                resi_in_chain.append(resi)
        # remove everything just before and after those residues
        max_resi = max(resi_in_chain)
        if max_resi != pose.size():
            pose.delete_residue_range_slow(max_resi + 1, pose.size())
        min_resi = min(resi_in_chain)
        if min_resi != 1:
            pose.delete_residue_range_slow(1, min_resi - 1)
        return pose

    def read_pickle(self, fin):
        """Read a pickle file. """
        with open(fin, 'rb') as f:
            data_pkl = pickle.load(f)
        return data_pkl

    def read_json(self, fin):
        """Read a JSON file."""
        with open(fin, 'r') as f:
            data_json = json.load(f)
        return data_json

    def ca_contact_matrix(self, pose):
        ca_xyz = []
        for resi in range(1, pose.size() + 1):
            ca_xyz.append(pose.residue(resi).atom("CA").xyz())
        ca_xyz = np.array(ca_xyz)
        return cdist(np.array(ca_xyz), np.array(ca_xyz))

    def find_disconnected_residues(self, pose, nieghbours_away = 10, distance_cutoff = 8):
        """Finds residues that are not in contact ('disconnected') with other residues a certain amount of residues upstream/downstream of it.

        :param nieghbours_away: The number of residues to consider in contact with a residue upstream and downstream of it in the sequence.
        :param distance_cutoff: The distance (in Angstrom) for which 2 atoms must be in contact in order for their respective residues
        to be considered in contact.
        :return a list of bools designating which residues are disconnected i.e. are not in contact (True) and in contact (False)."""
        cm, resi_map = contactmap_with_mapped_resi(pose, distance_cutoff=distance_cutoff)
        disconnected = np.full(pose.size(), False)
        for r0 in range(pose.size() - nieghbours_away):
            ai_start, ai_stop = resi_map[r0]["start"], resi_map[r0]["stop"]
            upstream, downstream = r0 + nieghbours_away, r0 - nieghbours_away
            upstream_contact, downstream_contact = False, False
            if upstream in resi_map:
                upstream_contact = np.any(cm[ai_start:ai_stop + 1, resi_map[upstream]["start"]:])
            if downstream in resi_map:
                downstream_contact = np.any(cm[ai_start:ai_stop + 1, :resi_map[downstream]["stop"] + 1])
            disconnected[r0] = not upstream_contact and not downstream_contact
        # get_pymol_selections(disconnected, pose)
        return disconnected

    def get_cut_points(self, pdb_info):
        """find cut points from N-termini and C-termini"""
        plddt = pdb_info["resi_plddt"]
        is_predicted_l = np.array(plddt).mean(axis=0) >= self.plddt_confidence
        dssp = [list(i) for i in pdb_info["dssp"]]
        is_l_l = ~(((np.char.equal(dssp, "H") + np.char.equal(dssp, "E")).sum(axis=0) / len(dssp) * 100) >= self.ss_confidence)
        disconnected = pdb_info["disconnected"]
        is_disconnected_l = np.array(disconnected).mean(axis=0) >= self.disconnect_confidence
        n_resi, c_resi = 0, 0
        n_term_loop = zip(is_predicted_l, is_l_l, is_disconnected_l)
        c_term_loop = zip(reversed(is_predicted_l), reversed(is_l_l), reversed(is_disconnected_l))
        for loop, cut_dir in zip((n_term_loop, c_term_loop), ("N", "C")):
            for is_predicted_well, is_l, is_disconnected in loop:
                # Old condition. This was really hard to read so I made it easier to understand
                # condition_1 = is_l and disconnected
                # condition_2 = is_l and not is_predicted_well
                # if not condition_1 and not condition_2:
                #     break
                connected = not disconnected
                good_prediction = is_predicted_well
                is_H_or_E_segment = not is_l
                if good_prediction and connected and is_H_or_E_segment:
                    break
                if cut_dir == "N":
                    n_resi += 1
                else:
                    c_resi += 1
        return n_resi, c_resi

    @staticmethod
    def cut_monomeric_pose(pose, n_resi, c_resi):
        """Cut a pose with a single chain"""
        if pose.size() == n_resi or pose.size() == c_resi:
            print("Will not cut pose as n_resi or c_resi is equal to the pose length")
            return pose
        # this cuts including the first and last index
        if c_resi > 0:
            pose.delete_residue_range_slow(pose.size() - c_resi + 1, pose.size())
        if n_resi > 0:  # else keep the N termini
            pose.delete_residue_range_slow(1, n_resi)
        return pose

    def construct_pdb_info(self):
        """Constructs an empty pdb_info dictionary object to be filled in during the apply call."""
        return {"model_name": [], "avg_plddt": [], "resi_plddt": [], "model_n": [], "ranking": [], "prediction_n": [], "mer": [],
                    "iptm+ptm":[], "ptm": [], "iptm": [], "unique_id": [],
                    "chain": [], "multimer":[], "dssp": [], "disconnected": [], "monomeric_poses": [], "monomeric_poses_cut": []}

    def sort_multimer_models(self, models):
        """sort the models so that they appear in the following order:
        model_1_*_0, model_2_*_0 ... then model_1_*_1, model_2_*_0 ... then model_1_*_n, model_2_*_n
        """
        # 1 first but each model type into its own container
        if self.model_type == "ranked":
            model_ns = {k:[] for k in map(str, range(0, len(models)))}
        else:
            model_ns = {k:[] for k in map(str, range(1, len(models) + 1))}
        for p in models:
            model_ns[p.name.split("_")[-1].split(".pdb")[0]].append(p)
        # 2 now sort each container based on the output number
        fs = lambda p: int(p.name.split("_")[-1].split(".")[0])
        model_ns = {k: sorted(v, key=fs) for k, v in model_ns.items()}
        out_models = []
        for vv in zip(*model_ns.values()):
            for v in vv:
                out_models.append(v)
        return out_models

    def get_model_types(self, prediction, orders):
        """Get the paths to the models to be used."""
        model_paths = list(prediction.glob(f"{self.model_type}*"))
        if len(model_paths) == 0:
            raise ValueError(f"No {self.model_type} models found in {prediction}!")
        return model_paths

    def get_multimer_predictions(self, pdbid, pdb_info):
        """Gets the multimer predictions for the pdbid and fills it in into the pdb_info object."""
        multimer_predictions = [p for p in Path(self.af_output_path).glob(f"{pdbid}_*_*") if not "_1_" in str(p.name)]
        if len(multimer_predictions) > 0:
            for completed, prediction in enumerate(multimer_predictions):
                if not has_all_files(prediction):
                    print(pdbid, "does not have all files!")
                    continue
                mers = int(prediction.name.split("_")[1])
                ranking_debug = self.read_json(prediction.joinpath("ranking_debug.json"))
                iptm_plus_ptms = ranking_debug["iptm+ptm"].values()
                orders = ranking_debug["order"]
                # find the number of predictions
                # num_models = len(list(prediction.glob(f"relaxed_model_1_multimer_v2_pred_*.pdb")))
                # models = [prediction.joinpath(f"relaxed_model_{i}_multimer_v2_pred_{j}.pdb") for i in range(1, 6) for j in range(num_models)]
                # get relaxed models if we have them, else take the unrelaxed
                model_paths = self.get_model_types(prediction, orders)
                models = self.sort_multimer_models(model_paths)
                # reduce the models if max_models is set
                models_accepted = 0
                for model, iptm_plus_ptm, order in zip(models, iptm_plus_ptms, orders):
                    resi_plddts = self.read_pickle(prediction.joinpath(f"result_{order}.pkl"))["plddt"]
                    resi_plddts_per_chain = np.array_split(resi_plddts, mers)
                    assert len(set([len(i) for i in resi_plddts_per_chain])) == 1
                    iptm = self.read_pickle(prediction.joinpath(f"result_{order}.pkl"))["iptm"]
                    ptm = self.read_pickle(prediction.joinpath(f"result_{order}.pkl"))["ptm"]
                    a_chain_was_accepted = False
                    for resi_plddt, chain in zip(resi_plddts_per_chain, ["A", "B", "C", "D", "E"]):
                        avg_plddt = resi_plddt.mean()
                        print("plddt:", avg_plddt, "ptmplus:", iptm_plus_ptm)
                        if (iptm_plus_ptm >= self.iptm_plus_ptm_threshold) and (avg_plddt >= self.avg_plddt_threshold):
                            a_chain_was_accepted = True
                            pose = self.cut_all_but_chain(pose_from_file(str(model)), chain)
                            pdb_info["monomeric_poses"].append(pose)
                            pdb_info["dssp"].append(Dssp(pose).get_dssp_secstruct())
                            pdb_info["disconnected"].append(self.find_disconnected_residues(pose).tolist())
                            pdb_info["avg_plddt"].append(avg_plddt)
                            pdb_info["model_name"].append(model)
                            pdb_info["model_n"].append(order.split("model_")[1].split("_")[0])
                            pdb_info["resi_plddt"].append(resi_plddt.tolist())
                            pdb_info["ranking"].append(orders.index(order) + 1)
                            # TODO: fix the options in this script. have a cn option that sets what the oligomer type it is we are trying to extract
                            if prediction.name.count("_") > 1:
                                mer = prediction.name.split("_")[-2]
                            else:
                                mer = prediction.name.split("_")[-1]
                            pdb_info["mer"].append(mer)
                            prediction_n = model.parent.name.split("_")[-1]
                            pdb_info["prediction_n"].append(prediction_n)
                            pdb_info["chain"].append(chain)
                            pdb_info["iptm+ptm"].append(iptm_plus_ptm)
                            pdb_info["iptm"].append(iptm)
                            pdb_info["ptm"].append(ptm)
                            pdb_info["multimer"].append(True)
                            pdb_info["unique_id"].append(f"{model.stem}_{chain}_{mer}_{prediction_n}_multi")
                    if a_chain_was_accepted:
                        models_accepted += 1
                    if self.max_multimers is not None and models_accepted == self.max_multimers:
                        return

    def sort_monomer_models(self, models):
        """sort the models so that they appear in the following order:
        <PDBID>_1_1/* <PDBID>_1_2/* <PDBID>_1_3 ...
        """
        return sorted(models, key=lambda p: int(p.name.split("_")[-1]))

    def get_monomer_predictions(self, pdbid, pdb_info):
        """Gets the monomer predictions for the pdbid and fills it in into the pdb_info object."""
        monomer_predictions = self.sort_monomer_models([p for p in Path(self.af_output_path).glob(f"{pdbid}_1_*")])
        # monomer predictions
        if len(monomer_predictions) > 0:
            for completed, prediction in enumerate(monomer_predictions):
                if not has_all_files(prediction):
                    continue
                ranking_debug = self.read_json(prediction.joinpath("ranking_debug.json"))
                orders = ranking_debug["order"]
                avg_plddts = ranking_debug["plddts"].values()
                models = [prediction.joinpath(f"relaxed_model_{i}_pred_0.pdb") for i in range(1, 6)]
                resi_plddts = [self.read_pickle(prediction.joinpath(f"result_model_{i}_pred_0.pkl"))["plddt"] for i in range(1, 6)]
                models_accepted = 0
                for n, (avg_plddt, model, resi_plddt) in enumerate(zip(avg_plddts, models, resi_plddts), 1):
                    # print("plddt:", avg_plddt)
                    if avg_plddt >= self.avg_plddt_threshold:
                        models_accepted += 1
                        pose = pose_from_file(str(model))
                        pdb_info["monomeric_poses"].append(pose)
                        pdb_info["dssp"].append(Dssp(pose).get_dssp_secstruct())
                        pdb_info["disconnected"].append(self.find_disconnected_residues(pose).tolist())
                        pdb_info["avg_plddt"].append(avg_plddt)
                        pdb_info["model_name"].append(model)
                        pdb_info["model_n"].append(n)
                        pdb_info["resi_plddt"].append(resi_plddt.tolist())
                        pdb_info["ranking"].append(orders.index(f"model_{n}_pred_0") + 1)
                        mer, prediction_n = prediction.name.split("_")[1:]
                        pdb_info["mer"].append(mer)
                        pdb_info["prediction_n"].append(prediction_n)
                        pdb_info["chain"].append("A")
                        pdb_info["iptm+ptm"].append(np.nan)
                        pdb_info["iptm"].append(np.nan)
                        pdb_info["ptm"].append(np.nan)
                        pdb_info["multimer"].append(False)
                        pdb_info["unique_id"].append(f"{model.stem}_{prediction_n}_{mer}_mono")
                        if self.max_monomers is not None and models_accepted == self.max_monomers:
                            return

    def reduce_set(self, pdb_info):
        """Reduces the information in the pdb_info object based on the pairwise RMSD of the cut_poses."""
        rmsd_diff = 0.005
        rmsd_threshold = self.init_rmsd_threshold + rmsd_diff
        if self.modify_rmsd_to_reach_min_models is None:
            attempts, minimum_models = 1, 1
        else:
            attempts = 18 # if start 0.1 it will end at 0.01 with -= 0.005 per iteration
            minimum_models = self.modify_rmsd_to_reach_min_models
        for attempt in range(attempts):
            rmsd_threshold -= rmsd_diff
            poses = pdb_info["monomeric_poses_cut"]
            # save matrix for later use
            rmsd_matrix = []
            for n, pose_i in enumerate(poses, 1):
                rmsd_ij = []
                for pose_j in poses:
                    rmsd_ij.append(CA_rmsd(pose_i, pose_j))
                rmsd_matrix.append(rmsd_ij)
                print(f"Calculated {n}/{len(poses)} RMSDs")
            indices_to_delete = []
            for i, (row, i_avg_plddt) in enumerate(zip(rmsd_matrix, pdb_info["avg_plddt"])):
                row = np.array(row)
                indices_under_threshold = [j for j in np.where((row < rmsd_threshold))[0] if j != i]
                for j in indices_under_threshold:
                    if pdb_info["avg_plddt"][j] < i_avg_plddt:
                        indices_to_delete.append(i)
                        break
            models_left = len(poses) - len(indices_to_delete)
            if models_left >= minimum_models:
                print(f"Reached the minimum number of models: {minimum_models}")
                break
        new_pdb_info = {k: [v for n, v in enumerate(vv) if n not in indices_to_delete] for k, vv in pdb_info.items()}
        assert len(new_pdb_info["monomeric_poses_cut"]) == models_left
        assert len(new_pdb_info["avg_plddt"]) == models_left
        print(f"{len(indices_to_delete)}/{len(poses)} poses were removed!")
        return new_pdb_info, rmsd_threshold

    def get_aligned_info(self, align_structure, pdb_info, n_resi, c_resi):
        """Store alignment information (TMscore and RMSD ) and if superimpose them"""
        # get the aligned structure as a pose object it and cut the same as the ensemble
        align_pose = pose_from_file(str(align_structure))
        align_pose_cut = self.cut_monomeric_pose(align_pose.clone(), n_resi, c_resi)
        # not align all the other structures to that. Also calculate the TMscore and RMSD.
        cut_str = f"{self.plddt_confidence}|{self.ss_confidence}|{self.disconnect_confidence}"
        pdb_info.update({"poses_aligned": [], "cut_poses_aligned": []})
        pdb_info.update({"tmscore": [], f"tmscore@{cut_str}": [], "rmsd": [], f"rmsd@{cut_str}": []})
        for pose_af, pose_af_cut in zip(pdb_info["monomeric_poses"], pdb_info["monomeric_poses_cut"]):
            # calculate TMscore / RMSD
            pdb_info["tmscore"].append(tmalign(pose_af.clone(), align_pose))
            pdb_info["rmsd"].append(CA_rmsd(align_pose, pose_af))
            pdb_info[f"tmscore@{cut_str}"].append(tmalign(pose_af_cut.clone(), align_pose_cut))
            pdb_info[f"rmsd@{cut_str}"].append(CA_rmsd(align_pose_cut, pose_af_cut))
            # superimpose poses
            pose_af_align = pose_af.clone()
            pose_af_cut_align = pose_af_cut.clone()
            calpha_superimpose_pose(pose_af_align, align_pose)
            calpha_superimpose_pose(pose_af_cut_align, align_pose_cut)
            pdb_info["poses_aligned"].append(pose_af_align)
            pdb_info["cut_poses_aligned"].append(pose_af_cut_align)

    @staticmethod
    def cut_multimer_poses(multimeric_pose, n_resi, c_resi):
        multimeric_pose_cut = None
        for pose in multimeric_pose.split_by_chain():
            pose_cut = AF2EvoDOCK.cut_monomeric_pose(pose, n_resi, c_resi)
            if multimeric_pose_cut is None:
                multimeric_pose_cut = pose_cut
            else:
                multimeric_pose_cut.append_pose_by_jump(pose_cut, multimeric_pose_cut.chain_end(multimeric_pose_cut.num_chains()))
        return multimeric_pose_cut

    def construct_if_not_symmetric(self, pose, chains_allowed, avg_plddt, cn):
        fold, _ = SymmetryMapper.get_fold_and_vrt_id(self.symmetry, cn)
        cs = CubicSetup()
        cs.load_norm_symdef(self.symmetry, fold)

        # get distance between the CA atoms
        com_subunits = []
        for i in range(1, pose.num_chains() + 1):
            com_subunits.append(pose.residue(residue_center_of_mass(pose, pose.chain_begin(i), pose.chain_end(i))).atom("CA").xyz())
        com_complex = np.array(center_of_mass(pose, 1, pose.size()))
        com_subunits = np.array(com_subunits)
        dst = []
        for com_subunit in com_subunits:
            dst.append(np.linalg.norm(com_subunit - com_complex))
        x_trans = np.array(dst).mean()

        # get the best chain
        best_chain = chains_allowed[np.argmax(avg_plddt)]

        # center pose
        best_pose = None
        for p in pose.split_by_chain():
            chain = p.pdb_info().chain(1)
            if chain == best_chain:
                best_pose = p
                break
        assert best_pose is not None
        best_pose.center()
        return cs, best_pose, x_trans, best_chain

        # flip it


    def globalfrommultimer(self, pdb_info, cn, n_resi, c_resi):
        # pdb_info_new = {k:[] for k in pdb_info.keys()}
        pdb_info_new = {"model_name": [], "x_trans": [], "chain": [], "input": [], "input_flip": []}
        models = set(pdb_info["model_name"])
        avg_plddts = [[c for c, m in zip(pdb_info["avg_plddt"], pdb_info["model_name"]) if model == m] for model in models]
        allowed_chains = [[c for c, m in zip(pdb_info["chain"], pdb_info["model_name"]) if model == m] for model in models]
        sm = SymmetryMapper()
        # align
        for model, chains_allowed, avg_plddt in zip(models, allowed_chains, avg_plddts):
            multimeric_pose_cut = self.cut_multimer_poses(pose_from_file(str(model)), n_resi, c_resi)
            cs, input_pose, input_pose_flip, input_pose_asym, input_pose_flip_asym, x_trans, chain_used = sm.run(model=multimeric_pose_cut, cn=cn,
                                                                                                                 symmetry=self.symmetry,
                                                                                                                 chains_allowed=chains_allowed,
                                                                                                                 T3F=False)  # fixme: remove T3F for public code
            if cs is None:
                if self.use_only_symmetric:
                    continue
                else:
                    # create symmetry
                    cs, input_pose_asym, x_trans, chain_used = self.construct_if_not_symmetric(multimeric_pose_cut, chains_allowed, avg_plddt, cn)
                    input_pose_flip_asym = input_pose_asym.clone()
                    sm.do_a_180_around_axis(input_pose_flip_asym, [1,0,0])
                    # cs, input_pose, input_pose_flip, input_pose_asym, input_pose_flip_asym, x_trans = \
            pdb_info_new["model_name"].append(model)
            pdb_info_new["x_trans"].append(x_trans)
            pdb_info_new["chain"].append(chain_used)
            if self.is_direction_allowed("up"):
                pdb_info_new["input"].append(input_pose_asym)
            if self.is_direction_allowed("down"):
                pdb_info_new["input_flip"].append(input_pose_flip_asym)

        if len(pdb_info_new["model_name"]) == 0:
            raise ValueError("No symmetrical models were found!")

        # Now we clean the pdb_info
        indices_to_delete = []
        for n, (model, chain) in enumerate(zip(pdb_info["model_name"], pdb_info["chain"]), 0):
            delete = True
            for model_new, chain_new in zip(pdb_info_new["model_name"], pdb_info_new["chain"]):
                # if model and chain matches to what already exists in pdb_info_new, then add it
                if model == model_new and chain == chain_new:
                    delete = False
                    break
            if delete:
                indices_to_delete.append(n)

        pdb_info_new.update(
            {k: [v for n, v in enumerate(vv) if n not in indices_to_delete] for k, vv in pdb_info.items() if k != "model_name"})

        return pdb_info_new

    def apply(self, model_name, align_structure=None, cn=None):
        # construct an emtpy pdb_info dictionary containing all information we want to store
        pdb_info = self.construct_pdb_info()

        # make a subdirectory in which you will dump the output files
        pdb_out = self.out_dir.joinpath(model_name).joinpath("pdbs")
        data_out = self.out_dir.joinpath(model_name).joinpath("data")
        pdb_out.mkdir(exist_ok=True, parents=True)
        data_out.mkdir(exist_ok=True, parents=True)

        # fill into the pdb_info object information from either or both the monomer or multimeric predictions
        if self.ensemble != "GlobalFromMultimer":
            self.get_monomer_predictions(model_name, pdb_info)
        self.get_multimer_predictions(model_name, pdb_info)
        assert len(set(pdb_info["unique_id"])) == len(pdb_info["model_name"]), "Not all names are unique - this should not happen!"

        # Checks if we have found accepted predictions. Returns an Error if not.
        if not pdb_info['model_name']:  # empty
            # --avg_plddt_threshold', help="The minimum value of the average pLDDT needed for the structure to be included in the"
            #                           "ensemble.", type=float, default=90)
            #     parser.add_argument('--iptm_plus_ptm_threshold'
            raise ValueError("No AlphaFold predictions have been accepted. Consider lowering --avg_plddt_threshold and --iptm_plus_ptm_threshold or run more AlphaFold predictions.")

        # Use the information in the pdb_info object obtain the N and C termini cut pointS
        n_resi, c_resi = self.get_cut_points(pdb_info)

        # Create another entry in the pdb_info object with all the poses cut
        pdb_info["monomeric_poses_cut"] = [self.cut_monomeric_pose(pose.clone(), n_resi, c_resi) for pose in pdb_info["monomeric_poses"]]

        # reduce the set of poses / other information based on RMSD in the pdb_info object.
        pdb_info, final_rmsd = self.reduce_set(pdb_info)

        # We need to symmetrize the models and retrieve data from them if the ensemble is GlobalFromMultimer
        if self.ensemble == "GlobalFromMultimer":
            pdb_info = self.globalfrommultimer(pdb_info, cn, n_resi, c_resi)

        # if we have more than self.total_max_models filter based on plddt
        self.filter_to_max(pdb_info)

        # store alignment information (TMscore and RMSD) and store aligned poses
        if align_structure is not None:
            self.get_aligned_info(align_structure, pdb_info, n_resi, c_resi)

        # prepack all the poses if the prepacking option is set
        score_data = None
        if self.prepack:
            score_data = self.pack_poses(pdb_info, n_resi, c_resi, align_structure)

        # output the final ensemble
        if self.ensemble == "GlobalFromMultimer":
            self.output_globalfrommultimer(pdb_info, pdb_out, model_name, cn, data_out, score_data)
        elif self.ensemble == "Local":
            self.output_local(pdb_info, pdb_out, data_out, model_name, score_data)

        self.output_data(pdb_info, data_out, model_name, n_resi, c_resi, final_rmsd)

    @staticmethod
    def cn_to_sympath(symmetry, cn):
        sympath = f"{symmetry}/{symmetry}_{{}}_norm.symm"
        data = Path(DATA)
        if symmetry == "I":
            if cn == "5":
                sympath = data.joinpath(sympath.format("HF"))
            else:
                sympath = data.joinpath(sympath.format(f"{cn}F"))
        elif symmetry == "O":
            if cn == "4":
                sympath = data.joinpath(sympath.format("HF"))
            else:
                sympath = data.joinpath(sympath.format(f"{cn}F"))
        elif symmetry == "T":
            if cn == "3":
                sympath = data.joinpath(sympath.format("HF"))
            else:
                sympath = data.joinpath(sympath.format(f"{cn}F"))
        return sympath


    def create_sample_model(self, pose_name, flip, cn, x_trans):
        sympath = self.cn_to_sympath(self.symmetry, cn)
        pose = pose_from_file(str(pose_name))
        init("-symmetry:initialize_rigid_body_dofs")
        SetupForSymmetryMover(str(sympath)).apply(pose)
        jid = CubicSetup.get_jumpidentifier_from_pose(pose)
        set_jumpdof_str_str(pose, f"JUMP{jid}fold111", "x", x_trans)
        from cubicsym.utilities import add_id_to_pose_w_base
        add_id_to_pose_w_base(pose, "1")
        ccsc = CloudContactScoreContainer(pose, str(sympath))
        icss = InitCubicSymmetrySlider(pose, str(sympath), ccsc)
        foldid = CubicSetup.get_jumpidentifier_from_pose(pose)
        icss.ccsc.set_ccs_and_cmc(pose)
        icss.slide_away(pose, foldid)
        pose.dump_pdb(str(pose_name.parent.parent.joinpath(f"{flip}_MODEL_of_{pose_name.name}")))

    def output_data(self, pdb_info, data_out, pdbid, n_resi, c_resi, final_rmsd):
        """Outputs data."""
        # first remove all pose info as we dont need them anymore
        for k in [i for i in pdb_info.keys() if "pose" in i or "input" in i]:
            del pdb_info[k]

        # output the pdb_info
        plddts = pdb_info.pop("resi_plddt")
        plddts_df = pd.DataFrame(columns=range(1, len(plddts[0]) + 1), data=np.array(plddts))
        pdb_df = pd.concat([pd.DataFrame(pdb_info), plddts_df], axis=1)
        pd.DataFrame(pdb_df).to_csv(data_out.joinpath(f"{pdbid}.csv"), index=False)


        # in meta also store how many monomer prediction vs multimer predictions got included
        pd.DataFrame(index=[0], data={"plddt_confidence": self.plddt_confidence,
                                      "ss_confidence": self.ss_confidence,
                                      "disconnect_confidence": self.disconnect_confidence,
                                      "n_resi": n_resi, "c_resi": c_resi,
                                      "avg_plddt_threshold": self.avg_plddt_threshold,
                                      "iptm_plus_ptm_threshold": self.iptm_plus_ptm_threshold,
                                      "direction": self.direction,
                                      "ensemble_rmsd_threshold": final_rmsd}).to_csv(data_out.joinpath(f"{pdbid}_meta.csv"), index=False)


    def filter_to_max(self, pdb_info):
        """Filters the pdb_info object so that the amount of elements per key does not exceed self.total_max_models. The structures are
        filtered based on the avg_plddt value so that the prediction with the best avg_plddt values are included."""
        avg_plddt = pdb_info["avg_plddt"]
        if self.total_max_models is not None and len(avg_plddt) > self.total_max_models:
            sorter = [i for i in reversed(np.argsort(avg_plddt))][:self.total_max_models]
            pdb_info_new = {k:[] for k in pdb_info.keys()}
            for k, vv in pdb_info.items():
                for i in sorter:
                    pdb_info_new[k].append(vv[i])
            return pdb_info_new
        else:
            return pdb_info

    def __pack_subroutine(self, pose, scores, fa_scfxn, packed_poses, align_structure=None):
        pose.conformation().detect_disulfides()
        scores["init_score"].append(fa_scfxn.score(pose))
        pack_pose = PosePacker(pose, "pdb", mode="custom").pack()
        scores["pack_score"].append(fa_scfxn.score(pack_pose))
        if align_structure is not None:
            scores["rmsd"].append(all_atom_rmsd(align_structure, pack_pose))
        packed_poses.append(pack_pose)

    def pack_poses(self, pdb_info, n_resi, c_resi, align_structure = None):
        """Packs all the relevant poses in the pdb_info object."""
        init(extra_options=opts) # initialize Rosetta as is done in prepacking.py
        scores = {"model":[], "init_score":[], "pack_score":[] } # model is filled in during output of the poses
        if align_structure is not None:
            scores["rmsd"] = []
            align_structure = self.cut_monomeric_pose(pose_from_file(str(align_structure)), n_resi, c_resi)
        fa_scfxn = get_fa_scorefxn()
        if self.ensemble == "GlobalFromMultimer":
            for flip, direction in zip(('up', 'down'), ("input", "input_flip")):
                if self.is_direction_allowed(flip):
                    packed_poses = []
                    for pose in pdb_info[direction]:
                        self.__pack_subroutine(pose, scores, fa_scfxn, packed_poses, align_structure)
                    pdb_info[direction] = packed_poses
        elif self.ensemble == "Local":
            packed_poses = []
            for pose in pdb_info["cut_poses_aligned"]:
                self.__pack_subroutine(pose, scores, fa_scfxn, packed_poses, align_structure)
            pdb_info["cut_poses_aligned"] = packed_poses
        return scores

    def output_local(self, pdb_info, pdb_out, data_out, pdbid, score_data=None):
        for pose, unique_id in zip(pdb_info["cut_poses_aligned"], pdb_info["unique_id"]):
            pose_name = str(pdb_out.joinpath(f"input_{unique_id}{'.prepack' if self.prepack else ''}.pdb"))
            pose.dump_pdb(pose_name)
            score_data["model"].append(pose_name)
        self.output_score_data(score_data, data_out, pdbid)

    def is_direction_allowed(self, direction):
        """Check if the direction is allowed"""
        if self.direction == "both" and (direction == "up" or direction == "down"):
            return True
        else:
            return self.direction == direction

    def output_globalfrommultimer(self, pdb_info, pdb_out, pdbid, cn, data_out, score_data=None):
        # output models
        x_trans_data = {"model": [], "x_trans": []}
        pose_names = []
        for flip, direction in zip(('up', 'down'), ("input", "input_flip")):
            if self.is_direction_allowed(flip):
                Path(pdb_out).joinpath(flip).mkdir(exist_ok=True)
                for unique_id, pose, x_trans in zip(pdb_info["unique_id"], pdb_info[direction], pdb_info["x_trans"]):
                    pose_name = Path(f"{pdb_out}/{flip}/input_{unique_id}_{flip}{'.prepack' if self.prepack else ''}.pdb")
                    pose_names.append(pose_name.name)
                    pose.dump_pdb(str(pose_name))
                    x_trans_data["model"].append(pose_name.name)
                    x_trans_data["x_trans"].append(x_trans)
                # output the last model
                self.create_sample_model(pose_name, flip, cn, x_trans)
        score_data["model"] = pose_names
        pd.DataFrame(x_trans_data).to_csv(f"{data_out}/{pdbid}_{cn}_xtrans.csv", index=False)
        self.output_score_data(score_data, data_out, pdbid)

    def output_score_data(self, score_data, data_out, pdbid):
        pd.DataFrame(score_data).to_csv(f"{data_out}/{pdbid}_pack_info.csv", index=False)

def has_all_files(prediction):
    files = [p.name for p in Path(prediction).glob(f"*")]
    if not "ranking_debug.json" in files:
        return False
    if not [f for f in fnmatch.filter(files, "*pkl") if f != "features.pkl"]:
        return False
    # if not any([not "pred" in s for s in fnmatch.filter(files, "result_model*.pkl")]):
    #     return False
    # _pred_0.pkl
    return True

def check_correct_format(path, ensemble):
    # the format should be <name of model>_<number of subunits>_*, * can be whatever but cannot contain an '_'.
    # check that all
    names = set()
    cns = []
    for p in Path(path).glob("*"):
        try:
            name, cn, _ = str(p.name).split("_")
        except ValueError:
            raise ValueError("The AF prediction folder must be of the format '<name of model>_<number of subunits>_*' and must not contain '_' other than the two '_' specified as it's used as the seperator.")
        names.add(name)
        cns.append(cn)
    assert len(names), "The name of the model must stay consistent across all folders. The following names are different: " \
                       f"{', '.join(names)}"
    assert all((i.isdigit() for i in cns)), "The number of subunits must all be an integer"
    cn = None
    if ensemble == "GlobalFromMultimer":
        assert len(set(cns)) == 1, "When using --ensemble=GlobalFromMultimer only one type of oligomeric type is supported"
        cn = cns[0]
    model_name = names.pop()
    return model_name, cn

def main(path, symmetry, out_dir, plddt_confidence, ss_confidence, disconnect_confidence, avg_plddt_threshold, iptm_plus_ptm_threshold, rmsd_threshold,
         ensemble, max_monomers, max_multimers, align_structure, total_max_models, modify_rmsd_to_reach_min_models,
         direction, model_type, use_only_symmetric):
    model_name, cn = check_correct_format(path, ensemble)
    init()
    out_dir.mkdir(exist_ok=True)
    af2evo = AF2EvoDOCK(symmetry, path, out_dir, plddt_confidence, ss_confidence, disconnect_confidence, avg_plddt_threshold, iptm_plus_ptm_threshold, rmsd_threshold,
                        ensemble, max_monomers, max_multimers, total_max_models, modify_rmsd_to_reach_min_models, direction, model_type, use_only_symmetric)
    # get all predicted pdb ids
    print(f"Analysing {model_name}")
    af2evo.apply(model_name, align_structure=align_structure, cn=cn)

if __name__ == "__main__":
    description = """Turns AlphaFold 2.2.2 predictions into an ensemble for EvoDOCK."""

    import argparse
    parser = argparse.ArgumentParser()
    # parsing general options
    parser.add_argument('--path', help="Path to the predictions. It should point to a directory which contains output folders of AlphaFold2 and/or AlphaFold-Multimer predictions."
                                       "It expects a particular name format for the folders. They should look like the following:'<name of model>_<number of subunits>_*' where * can be anything. It must not contain "
                                       "'_' other than the two '_' specified as it's used as the seperator. For example 1STM_1_ or 1STM_5_whatever. The _X_ is used to identify"
                                       "what the oligomeric state is. When running --ensemble=GlobalFromMultimer all of the folders must be of the same oligomeric type.", type=str, required=True)
    parser.add_argument('--symmetry', help="Symmetry to model", type=str, choices=("T", "O", "I"), required=True)
    parser.add_argument('--out_dir', help="Output directory in which to store the output ensemble and the processed data directory. If the directory does not exists, it will be created.", type=str, required=True)
    parser.add_argument('--model_type', help="Which output models from Alphafold to use.", type=str, default="ranked", choices=["ranked", "relaxed", "unrelaxed"])
    # parsing options
    parser.add_argument('--ensemble', help="Ensemble type to create. "
                          "Local: Creates an ensemble used for Local assembly docking with EvoDOCK. The whole ensemble will be aligned ontoi "
                          "a structure specified with --align_structure. The --align_structure option must be set with this option. "
                          "GlobalFromMultimer: Creates an ensemble used for Global assembly docking with EvoDOCK based on AlphaFold Multimer predicitions. ",
                        type=str, choices=["Local", "GlobalFromMultimer"], required=True)
    parser.add_argument('--align_structure', help="When running --ensemble=Local this has point to the structure one wants to align to. It should be the input " \
                          "that has been produced by the symmetrization script: 'cubic_to_rosetta.py'. See https://github.com/Andre-lab/cubicsym.", type=str)
    parser.add_argument('--max_monomers', help="The maximum number of monomer prediction models to include in the ensemble.", type=int)
    parser.add_argument('--max_multimers', help="The maximum number of multimer prediction models to include in the ensemble. Notice that, "
                          "since a multimer model consists of multiple chains, more than one chain can be included in "
                          "the ensemble from the same prediction model, without counting as more than one towards the maximum.", type=int)
    parser.add_argument('--max_total_models', help="The maximum number of models (monomer and multimer) to include in the final ensemble. This will "
                           "be the final maximum output ensemble size. However, this will not nescesarily be the final "
                           "ensemble size because this depends on the --rmsd_threshold value and the input models "
                           "pLDDT and iptm+ptm values. If more models are accepted than allowed with --max_total_models, "
                           "the script will filter the models based on their avg_plddt and output the ones with the "
                           "best avg_plddt.", type=int)
    parser.add_argument('--modify_rmsd_to_reach_min_models', help="Will iteratively lower the rmsd_threshold until the minimum number of models "
                           " (given by the value passed to the flag) are included in the ensemble. However, this will "
                           " not necessarily be the final ensemble size because this depends on input models "
                           " pLDDT and iptm+ptm values.", type=int)
    parser.add_argument("--direction", help="Chose which directions of the prediction to output. Options are 'up', 'down' and 'both'. The "
                           " former 2 options outputs either the up or down predictions and both outputs both of them. ",
                           choices=('up', 'down', 'both'), default="both")
    # threshold options for ensemble inclusion
    parser.add_argument('--avg_plddt_threshold', help="The minimum value of the average pLDDT needed for the structure to be included in the "
                          "ensemble.", type=float, default=90)
    parser.add_argument('--iptm_plus_ptm_threshold', help="The minimum value of average iptm+ptm needed for the multimeric structure to be "
                          "included in the ensemble.", type=float, default=0.90)
    parser.add_argument('--rmsd_threshold', help="Minimum pairwise RMSD between all structures in the ensemble.", type=float, default=0.1)
    # cutting options
    parser.add_argument('--plddt_confidence', help="The threshold for the mean score of plddt that stops cutting the terminis "
                          "if a residue i designated 'L' by DSSP", type=float, default=90)
    parser.add_argument('--ss_confidence', help="The percentage of predictions that must NOT be designated 'L' in order to stup cutting the "
                          "terminis", type=float, default=70)
    parser.add_argument('--disconnect_confidence', help="The percentage of predictions that must be designated 'disconnected' in order to "
                          "continue the cutting if a residue is designated 'L'", type=float, default=70)
    parser.add_argument('--use_only_symmetric', help="If --ensemble='GlobalFromMultimer', will only use AFM predictions that Rosetta predicts "
                                                     "to be symmetric. In that case the values for the symmetric paramters are "
                                                     "set based on the predictions. If set to to False, a generic symmetry will be used and "
                                                     "the symmetric parameters will not be based on the predictions. Setting it to False "
                                                     "can be useful if the AFM interface predictions are bad.", type=strtobool, default=True)
    # plddt_confidence, ss_confidence, disconnect_confidence
    args = parser.parse_args()

    if args.ensemble.upper() == "LOCAL" and args.align_structure is None:
        raise ValueError("when running --ensemble=Local, --align_structure must be set.")

    main(args.path, args.symmetry, Path(args.out_dir), args.plddt_confidence, args.ss_confidence, args.disconnect_confidence,
         args.avg_plddt_threshold, args.iptm_plus_ptm_threshold, args.rmsd_threshold, args.ensemble,
         args.max_monomers, args.max_multimers, args.align_structure, args.max_total_models, args.modify_rmsd_to_reach_min_models,
         args.direction, args.model_type, args.use_only_symmetric)
