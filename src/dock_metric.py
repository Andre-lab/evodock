#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Metric Calculator class
@Author: Mads Jeppesen
@Date: 8/1/22
"""

from pyrosetta.rosetta.protocols.symmetric_docking import SymDockProtocol
import numpy as np
from pyrosetta.rosetta.numeric import xyzVector_double_t
from pyrosetta.rosetta.protocols.docking import calc_interaction_energy, calc_Irmsd
from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreFunctionFactory, rms_at_all_corresponding_atoms
from pyrosetta import Vector1
from pyrosetta.rosetta.protocols.scoring import Interface
from pyrosetta import AtomID
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID


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
        return CA_rmsd(self.native, pose)
        # return self.symdock.calc_irms(pose) - CA_rmsd_symmetric is called under the hood anyways

    def i_rmsd(self, pose):
        """Calculate interface RMSD."""
        return calc_Irmsd(self.native, pose)

    def interaction_energy(self, pose):
        """Calculate interface energy across jump 1 of the pose."""
        return calc_interaction_energy(pose, self.score_func, Vector1([1]))


class SymmetryDockMetric(DockMetric):
    """Calculates metrics related to docking such as RMSD, interface RMSD and interface score."""

    def __init__(self, native, score_func=None, bound=True, jump_ids: list = None, dof_ids: list = None, trans_mags: list = None):
        """Initializes a SymmetryDockMetric object.

        :param native: symmetric native pose.
        :param score_func: score function to use when calculating the interaction energy.
        :param bound: If true the RMSD of the pose to the native is always against the same chains. If not all chains will be compared
            to all chains in order to get the correct RMSD. This can be very expensive.
        :param jump_ids: jump ids of the jumps to move when calculating the interaction energy. Should match dof_ids and trans_mag.
        :param dof_ids: dof ids of the dofs to move when calculating the interaction energy. should match jump_ids and trans_mag.
        :param trans_mags: translation magnitudes to employ when calculating the interaction energy. should match jump_ids and dof_ids.
        """
        super().__init__(native, score_func)
        self.native = native
        self.bound = bound
        # member variables have to be initialized differently depending on bound = False or True,
        # as all the metrics are calculated differently in each case.
        if not bound:
            self.symdock = SymDockProtocol()
            self.symdock.set_native_pose(self.native)
            self.atommap = None
        else:
            self.symdock = None
            self.atommap = self.__set_inteface_atoms()
        if jump_ids or dof_ids or trans_mags:
            assert len(jump_ids) == len(dof_ids) and len(jump_ids) == len(trans_mags), "jump_ids, dof_ids and trans_mag must be of equal length"
        self.jump_ids = jump_ids
        self.dof_ids = dof_ids
        self.trans_mags = trans_mags

    def ca_rmsd(self, pose):
        """Calculate symmetric RMSD."""
        if self.bound:
            return CA_rmsd(self.native, pose)
        else:
            # CA_rmsd_symmetric(self.native, pose) // CA_rmsd_symmetric is called under the hood anaways
            return self.symdock.calc_rms(pose)

    def i_rmsd(self, pose):
        """Calculate symmetric interface RMSD."""
        if self.bound:
            return rms_at_all_corresponding_atoms(self.native, pose, self.atommap)
        else:
            return self.symdock.calc_Irms(pose)

    def interaction_energy(self, pose):
        """Calculate interface energy."""
        assert self.jump_ids and self.dof_ids and self.trans_mags, \
            "jump_ids, dof_ids and trans_mags must be defined in the constructor to use this method!"
        self.__apply_trans(pose, 1)
        monomeric_score = self.score_func.score(pose)
        self.__apply_trans(pose, -1)
        multimeric_score = self.score_func.score(pose)
        return multimeric_score - monomeric_score

    def __apply_trans(self, pose, slidedir):
        """Applies all set translational degrees of freedom."""
        for jumpid, dofid, transmag in zip(self.jump_ids, self.dof_ids, self.trans_mags):
            flexible_jump = pose.jump(jumpid)
            trans = np.asarray(flexible_jump.rt().get_translation())
            trans[dofid - 1] += transmag * slidedir
            flexible_jump.set_translation(xyzVector_double_t(*trans))
            pose.set_jump(jumpid, flexible_jump)
        pose.conformation().detect_disulfides() # else disulfide energy can blow up!

    def __set_inteface_atoms(self):
        """Sets which atoms are determined to be in the interface in the native structure."""
        # The Interface class uses the energy to determine interface residues so we have to score it first.
        self.score_func(self.native)
        interface = Interface()
        interface.distance(10.0) # 10.0 as in  protocols::docking::calc_Irmsd
        interface.calculate(self.native)

        # In pyrosetta FArrays dont work well so we cannot use rmsd_with_super_subset as is used in calc_Irmsd.
        # Instead we use rms_at_all_corresponding_atoms which under the hood calls the same functions as CA_rmsd does which in turn
        # should accomplish the same as rmsd_with_super_subset is trying to accomplish.
        # We have to construct a std::map< core::id::AtomID, core::id::AtomID >. We will only use Heavy atoms because
        # this is what is done in protocols::docking::calc_Irmsd.
        atommap = map_core_id_AtomID_core_id_AtomID()
        for ri in range(1, self.native.size() + 1):
            if interface.is_interface(ri):
                for ai in range(1, self.native.residue(ri).natoms() + 1):
                    if self.native.residue(ri).atom_type(ai).is_heavyatom():
                        atommap[AtomID(ai, ri)] = AtomID(ai, ri)
        return atommap

    # def get_atommap_as_pymol_str(self):
    #     string = "color red, "
    #     for atomid in self.atommap:
    #         ri, ai = atomid.rsd(), atomid.atomno()
    #         string += f"resi {ri} "

    def get_interface_residues_as_pymol_str(self, split_lines:int=None):
        ris = []
        for atomid in self.atommap:
            ri, ai = atomid.rsd(), atomid.atomno()
            pdb_r, pdb_c = self.native.pdb_info().pose2pdb(ri).split()
            ris.append(f"(chain {pdb_c} and resi {pdb_r})")
        if split_lines:
            return ["color red, " + " OR ".join(l) for l in np.array_split(ris, split_lines)]
        else:
            return "color red, " + " OR ".join(ris)



