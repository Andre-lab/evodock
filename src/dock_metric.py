#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Metric Calculator class
@Author: Mads Jeppesen
@Date: 8/1/22
"""
import random

from pyrosetta.rosetta.protocols.symmetric_docking import SymDockProtocol
import numpy as np
from pyrosetta.rosetta.numeric import xyzVector_double_t
from pyrosetta.rosetta.protocols.docking import calc_interaction_energy, calc_Irmsd
from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreFunctionFactory, rms_at_all_corresponding_atoms, CA_rmsd_symmetric
from pyrosetta import Vector1
from pyrosetta.rosetta.protocols.scoring import Interface
from pyrosetta import AtomID
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID
from cubicsym.alignment import sequence_alignment_on_chain_set
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from cubicsym.cubicsetup import CubicSetup
from pyrosetta.rosetta.protocols.scoring import Interface
from cubicsym.utilities import get_all_ca_atoms_slow
from scipy.spatial.distance import pdist
from pyrosetta.rosetta.core.pose import residue_center_of_mass
from copy import deepcopy
from warnings import warn
from symmetryhandler.reference_kinematics import get_dofs
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_jump_num
from pyrosetta.rosetta.protocols.symmetric_docking import SymDockProtocol
from pyrosetta.rosetta.utility import vector1_std_pair_unsigned_long_unsigned_long_t

class DockMetric:
    """Calculates metrics related to docking such as RMSD, interface RMSD and interface score."""

    def __init__(self, native, score_func=None):
        """Initializes a DockMetric object.

        :param native: symmetric native pose.
        :param score_func: score function to use when calculating the interaction energy.
        """
        self.native = native
        self.score_func = score_func if score_func else ScoreFunctionFactory.create_score_function("ref2015")

    def ca_rmsd(self, pose):
        """Calculate RMSD."""
        if self.native is not None:
            return CA_rmsd(self.native, pose)
        else:
            return -1
        # return self.symdock.calc_irms(pose) - CA_rmsd_symmetric is called under the hood anyways

    def interface_rmsd(self, pose):
        """Calculate interface RMSD."""
        if self.native is not None:
            return calc_Irmsd(self.native, pose)
        else:
            return -1

    def interaction_energy(self, pose):
        """Calculate interface energy across jump 1 of the pose."""
        return calc_interaction_energy(pose, self.score_func, Vector1([1]))


class SymmetricDockMetric(DockMetric):
    """Calculates metrics related to docking such as RMSD, interface RMSD and interface score."""

    def __init__(self, native, score_func=None):
        """Initializes a DockMetric object.

        :param native: symmetric native pose.
        :param score_func: score function to use when calculating the interaction energy.
        """
        super().__init__(native, score_func)
        self.trans_jump = None
        self.symdockproto = SymDockProtocol()
        self.symdockproto.set_native_pose(native)

    def get_translational_jump(self, pose):
        """Gets the translational dof. Usually """
        if self.trans_jump is None:
            for k, v in pose.conformation().Symmetry_Info().get_dofs().items():
                if v.allow_dof(1) or v.allow_dof(2) or v.allow_dof(3):
                    self.trans_jump = Vector1([k])
                    return self.trans_jump
            raise ValueError("No translational DOF found!")
        else:
            return self.trans_jump


    # CA_rmsd_symmetric
    def ca_rmsd(self, pose):
        """Calculate RMSD."""
        if self.native is not None:
            # old approach: CA_rmsd_symmetric(self.native, pose)
            return self.symdockproto.calc_rms(pose)
        else:
            return -1

    def interface_rmsd(self, pose):
        """Calculate interface RMSD."""
        if self.native is not None:
            # old aproach
            # self.trans_jump = self.get_translational_jump(pose)
            # pose.conformation().Symmetry_Info().get_dofs()
            # calc_Irmsd(pose, self.native, self.score_func, self.get_translational_jump(pose)) - does not work with symmetry
            raise NotImplementedError("This does not work for some reason!!!!")
            return self.symdockproto.calc_Irms(pose)
        else:
            return -1

    def interaction_energy(self, pose):
        """Calculate interface energy across jump 1 of the pose."""
        return self.symdockproto.calc_interaction_energy(pose) #calc_interaction_energy(pose, self.score_func, self.get_translational_jump(pose))

class CubicDockMetric:
    """Calculates metrics related to docking such as RMSD, interface RMSD and interface score."""

    def __init__(self, crystallic_native, input_pose, native_symdef, input_symdef, score_func=None,
                 jump_ids: list = None, dof_ids: list = None, trans_mags: list = None, use_map=None):
        """Initializes a SymmetryDockMetric object.

        :param input_pose:
        :param crystallic_native: symmetric native pose.
        :param score_func: score function to use when calculating the interaction energy.
        :param bound: If true the RMSD of the pose to the native is always against the same chains. If not all chains will be compared
            to all chains in order to get the correct RMSD. This can be very expensive.
        :param jump_ids: jump ids of the jumps to move when calculating the interaction energy. Should match dof_ids and trans_mag.
        :param dof_ids: dof ids of the dofs to move when calculating the interaction energy. should match jump_ids and trans_mag.
        :param trans_mags: translation magnitudes to employ when calculating the interaction energy. should match jump_ids and dof_ids.
        """
        self.crystallic_native = crystallic_native
        self.cubicsetup_native = None
        self.same_handedness = None
        self.cubicsetup_input = CubicSetup(input_symdef)
        # This is important if using a native
        if crystallic_native is not None:
            assert not is_symmetric(self.crystallic_native)
            self.cubicsetup_native = CubicSetup(native_symdef)
            self.same_handedness = self.cubicsetup_input.righthanded == self.cubicsetup_native.righthanded
            self.CA_atom_map = self.cubicsetup_input.construct_atom_map_any2hf(input_pose, self.crystallic_native,
                                                                               same_handedness=self.same_handedness,
                                                                               interface=False,
                                                                               predicate="ca", use_map=use_map)
            self.CA_interface_atom_map = self.cubicsetup_input.construct_atom_map_any2hf(input_pose, self.crystallic_native,
                                                                                         same_handedness=self.same_handedness,
                                                                                         interface=True,
                                                                                         predicate="ca", use_map=use_map)
        # construct CA chain map to be used for CA_rmsd calculations
        if jump_ids or dof_ids or trans_mags:
            assert len(jump_ids) == len(dof_ids) and len(jump_ids) == len(trans_mags), "jump_ids, dof_ids and trans_mag must be of equal length"
        self.jump_ids = jump_ids
        self.dof_ids = dof_ids
        self.trans_mags = trans_mags
        if score_func is None:
            self.score_func = ScoreFunctionFactory.create_score_function("ref2015")
        else:
            self.score_func = score_func
        self.max_com_distance = self.get_max_distance_between_ca_coms(input_pose)
        self.pose_com_ca_atomids = self.get_pose_com_atom_ids(input_pose)

    def get_max_distance_between_ca_coms(self, pose, touching_buffer=10):
        """Gets the maximum distance 2 chains of the pose must have between their center of mass (COM) CA atoms in order for
        the chains not to touch each other. This is calculated as the largest distance from a single chains COM ca atom
        to any other CA atom within that chain and then multiplied by 2. A buffer is (touching_buffer) is added to this
        distance"""
        all_ca_xyzs = get_all_ca_atoms_slow(pose)
        com_resi = residue_center_of_mass(pose, pose.chain_begin(1), pose.chain_end(1))
        com_ca_xyz = np.array(pose.xyz(AtomID(2, com_resi)))  # 2 = CA
        return max(np.linalg.norm(com_ca_xyz - i) for i in all_ca_xyzs) * 2 + touching_buffer

    def get_pose_com_atom_ids(self, pose):
        atomids = []
        for chain in range(1, pose.num_chains()): # not +1 because we dont want VRT
            resi = residue_center_of_mass(pose, pose.chain_begin(chain), pose.chain_end(chain))
            assert pose.residue(resi).name() != "VRT"
            atomids.append(AtomID(2, resi)) # 2 = CA
        return atomids

    def get_min_com_ca_distances(self, pose):
        xyzs = []
        for atomid in self.pose_com_ca_atomids:
            xyzs.append(pose.xyz(atomid))
        return min(pdist(xyzs))

    def ca_rmsd(self, pose):
        """Calculate cubic symmetric RMSD."""
        if self.crystallic_native is not None:
            return self.cubicsetup_input.rmsd_hf_map_with_atom_map(pose, self.crystallic_native, self.CA_atom_map)
        else:
            return -1

    def interface_rmsd(self, pose):
        """Calculate symmetric interface RMSD."""
        if self.crystallic_native is not None:
            return self.cubicsetup_input.rmsd_hf_map_with_atom_map(pose, self.crystallic_native, self.CA_interface_atom_map)
        else:
            return -1

    def interaction_energy(self, pose):
        """Calculates interface energy."""
        assert self.jump_ids and self.dof_ids and self.trans_mags, \
            "jump_ids, dof_ids and trans_mags must be defined in the constructor to use this method!"
        # We calculate the the multimeric score here and clone the pose now.
        # I've experienced this error when translating away and back and call self.score_func.score(pose):
        # File: ... src/core/scoring/etable/count_pair/CountPairFactory.cc:220
        # [ ERROR ] UtilityExitException
        # ERROR: Assertion `res2.is_bonded(res1)` failed.
        # And I can't seem to resolve it. I believe it has to do with disulfide bonds but even with disulfide bond
        # detection through pose.conformation().detect_disulfides() it still happens. Therefore I work on a
        # clone pose. Test cases are 1JH5 and 1X36.
        pose.conformation().detect_disulfides()
        multimeric_score = self.score_func.score(pose)
        pose_clone = pose.clone()
        init_trans_mags = deepcopy(self.trans_mags)
        trans_mags_used, no_touching_achieved = self._translate_away_until_no_touching(pose_clone, init_trans_mags)
        if no_touching_achieved:
            monomeric_score = self.score_func.score(pose_clone)
            # self._apply_trans(pose, -1, trans_mags_used)
            return multimeric_score - monomeric_score
        else:
            return 10**6

    def _translate_away_until_no_touching(self, pose, trans_mags, attempts=10):
        """translate the pose away until none of the chains are touching. Try 'attempts' times."""
        for attempt in range(1, attempts + 1):
            # translate away
            self._apply_trans(pose, 1, trans_mags)
            # check if the chains touch each other
            if self.get_min_com_ca_distances(pose) <= self.max_com_distance:
                warn(f"The chain of of the pose might be touching each other when sliding away at with the following dofs: "
                     f"{get_dofs(pose)}"
                     f"Attempts to modify the applied translations attempt={attempt}")
                # translate back
                self._apply_trans(pose, -1, trans_mags)
                # modify the trans_mag
                rnd_index = random.randint(0, len(trans_mags) - 1)
                trans_mags[rnd_index] += random.randint(1, 11) * self.max_com_distance
            else:
                return trans_mags, True
        # if all attempts have been attempted, slide away again. We should always leave in a slide away state
        # as outside of this function we will slide onto again.
        self._apply_trans(pose, 1, trans_mags)
        return trans_mags, False

    def _apply_trans(self, pose, slidedir, trans_mags):
        """Applies all set translational degrees of freedom."""
        for jumpid, dofid, transmag in zip(self.jump_ids, self.dof_ids, trans_mags):
            flexible_jump = pose.jump(jumpid)
            trans = np.asarray(flexible_jump.rt().get_translation())
            trans[dofid - 1] += transmag * slidedir
            flexible_jump.set_translation(xyzVector_double_t(*trans))
            pose.set_jump(jumpid, flexible_jump)
        # fixme: THIS NEEDS TO BE HERE!
        # pose.conformation().detect_disulfides(vector1_std_pair_unsigned_long_unsigned_long_t())
        pose.conformation().detect_disulfides()

    # def __set_inteface_atoms(self, input_pose):
    #     """Sets which atoms are determined to be in the interface in the native structure."""
    #     # The Interface class uses the energy to determine interface residues so we have to score it first.
    #     self.score_func(self.crystallic_native)
    #     interface = Interface()
    #     interface.distance(10.0) # 10.0 as in  protocols::docking::calc_Irmsd
    #     interface.calculate(self.crystallic_native)
    #
    #     if is_symmetric(input_pose):
    #         # VRT's have their own chain num so therefore -1
    #         assert self.crystallic_native.num_chains() == input_pose.num_chains() - 1, "chains between the input pose and native does not match"
    #     else:
    #         assert self.crystallic_native.num_chains() == input_pose.num_chains(), "chains between the input pose and native does not match"
    #
    #     # In pyrosetta FArrays dont work well so we cannot use rmsd_with_super_subset as is used in calc_Irmsd.
    #     # Instead we use rms_at_all_corresponding_atoms which under the hood calls the same functions as CA_rmsd does which in turn
    #     # should accomplish the same as rmsd_with_super_subset is trying to accomplish.
    #     # We have to construct a std::map< core::id::AtomID, core::id::AtomID >. We will only use Heavy atoms because
    #     # this is what is done in protocols::docking::calc_Irmsd.
    #     atommap = map_core_id_AtomID_core_id_AtomID()
    #     for chain in range(1, self.crystallic_native.num_chains() + 1):
    #         alignment = sequence_alignment_on_chain_set(self.crystallic_native, input_pose, (chain,), (chain,))
    #         native_resi = alignment[f"1_{chain}"]
    #         input_pose_resi = alignment[f"2_{chain}"]
    #         for ri_nat, ri_inp in zip(native_resi, input_pose_resi):
    #             if interface.is_interface(ri_nat):
    #                 for ai in range(1, self.crystallic_native.residue(ri_nat).natoms() + 1):
    #                     if self.crystallic_native.residue(ri_nat).atom_type(ai).is_heavyatom():
    #                         atommap[AtomID(ai, ri_nat)] = AtomID(ai, ri_inp)
    #     return atommap

    # def get_atommap_as_pymol_str(self):
    #     string = "color red, "
    #     for atomid in self.atommap:
    #         ri, ai = atomid.rsd(), atomid.atomno()
    #         string += f"resi {ri} "

    # def get_interface_residues_as_pymol_str(self, split_lines:int=None):
    #     ris = []
    #     for atomid in self.rmsd_interface_atommap:
    #         ri, ai = atomid.rsd(), atomid.atomno()
    #         pdb_r, pdb_c = self.crystallic_native.pdb_info().pose2pdb(ri).split()
    #         ris.append(f"(chain {pdb_c} and resi {pdb_r})")
    #     if split_lines:
    #         return ["color red, " + " OR ".join(l) for l in np.array_split(ris, split_lines)]
    #     else:
    #         return "color red, " + " OR ".join(ris)



