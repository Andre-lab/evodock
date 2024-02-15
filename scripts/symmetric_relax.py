#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run relax with symmetry and allows one to output the final symmetry file and the input file centered
@Author: Mads Jeppesen
@Date: 6/15/22
"""
import argparse
from pyrosetta.rosetta.utility import vector1_std_pair_unsigned_long_unsigned_long_t
from pyrosetta import init, pose_from_file, standard_task_factory, Pose
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from symmetryhandler.symmetrysetup import SymmetrySetup
from pathlib import Path
import pandas as pd
from src.dock_metric import CubicDockMetric, SymmetricDockMetric
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_jump_num
from cubicsym.cubicsetup import CubicSetup
import numpy as np
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent,
    InitializeFromCommandline,
    NoRepackDisulfides,
    RestrictToRepacking,
)
from pyrosetta import get_fa_scorefxn
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from argparse import RawTextHelpFormatter


class PosePacker:
    def __init__(self, pose, filename, mode="standard"):
        self.pose = pose
        if "pdb" in filename:
            self.output_name = filename.replace(".pdb", ".prepack.pdb")
        elif "cif" in filename:
            self.output_name = filename.replace(".cif", ".prepack.cif")
        else:
            raise ValueError(f"{filename} not understood. Only '.pdb' or '.cif' is understood")
        self.mode = mode
        self.scorefxn = get_fa_scorefxn()

class PosePackRotamers(PosePacker):
    def pack(self):
        pack_pose = Pose()
        pack_pose.assign(self.pose)

        local_tf = TaskFactory()
        local_tf.push_back(InitializeFromCommandline())
        # local_tf.push_back(IncludeCurrent())
        local_tf.push_back(RestrictToRepacking())
        local_tf.push_back(NoRepackDisulfides())
        conformer_full_repack = PackRotamersMover(self.scorefxn)
        conformer_full_repack.task_factory(local_tf)
        conformer_full_repack.apply(pack_pose)
        return pack_pose



# FIXME: delete this
dof_map = {
    "x" : 1,
    "y" : 2,
    "z" : 3,
    "angle_x" : 4,
    "angle_y" : 5,
    "angle_z" : 6,
    # FIXME: DELETE THIS IN FUTURE WHEN FIXED IN THE RIGIDBODYDOFADAPTIVEMOVER
    "x_angle": 4,
    "y_angle": 5,
    "z_angle": 6,
    ################
}

def calculate_metrics(pose, native, input_pose, native_symdef, input_symdef, use_rmsd_map=None):
    # get the translational dfo
    info = dict(pose.scores.items())  # all score terms
    if CubicSetup(input_symdef).is_cubic():
        jumpints, dofints, trans_mags = [], [], []
        foldid =CubicSetup.get_jumpidentifier_from_pose(pose)
        for jump, dof, trans in zip([f"JUMP{foldid}fold1", f"JUMP{foldid}fold111"], ["z", "x"], [2000, 1000]):
            jumpints.append(sym_dof_jump_num(pose, jump))
            dofints.append(dof_map[dof])
            trans_mags.append(trans)
        dockmetric = CubicDockMetric(native, input_pose, native_symdef=native_symdef, input_symdef=input_symdef,
                                     jump_ids=jumpints, dof_ids=dofints, trans_mags=trans_mags, use_map=use_rmsd_map)
    else:
        dockmetric = SymmetricDockMetric(native)
    info["energy"] = info.pop("total_score")
    info["rmsd"] = dockmetric.ca_rmsd(pose)
    # fixme: interaction energy only works for cubic symmetry
    if CubicSetup(input_symdef).is_cubic():
        info["Irmsd"] = dockmetric.interface_rmsd(pose)
    else:
        info["Irmsd"] = None
    info["Ienergy"] = dockmetric.interaction_energy(pose)
    return info

def output_info(pose, native, input_pose, info_out, native_symdef, input_symdef, relaxes_done, use_old_data = None, use_rmsd_map=None):
    # get the translational dfos
    info = calculate_metrics(pose, native, input_pose, native_symdef, input_symdef, use_rmsd_map)
    if use_old_data is not None:
        # owverwrite
        info = {k: np.NaN for k, v in info.items()}
        info["rmsd"] = use_old_data["rmsd"].values[0]
        info["Irmsd"] = use_old_data["Irmsd"].values[0]
        info["energy"] = use_old_data["score"].values[0]
        info["Ienergy"] = use_old_data["Iscore"].values[0]
    for k, v in relaxes_done.items():
        info[k] = v
    pd.DataFrame(info, index=[0]).to_csv(Path(info_out), index=False)

def make_fastrelax(pose_file, cycles=5, jumps=True, bb=True, sc=True, constrain_coords=False, ex1=True, ex2aro=True, ex2=True):
    # we initialize Rosetta here and this sets some values inside Rosetta relax that we cannot set otherwise
    # (as far as I have looked!. At least this is the easiest)
    init_rosetta(pose_file, constrain_coords=constrain_coords, ex1=ex1, ex2aro=ex2aro, ex2=ex2)
    # initialize fastrelax with ref2015 and cycles
    sfxn = ScoreFunctionFactory.create_score_function("ref2015")
    fastrelax = FastRelax(sfxn, cycles)

    # Set the packer task
    # This is an exact version of the packer task inside the Fastrelax code
    tf = standard_task_factory()
    tf.push_back(RestrictToRepacking())
    tf.push_back(IncludeCurrent())
    tf.push_back(InitializeFromCommandline())
    fastrelax.set_task_factory(tf)

    # Set the movemap
    mp = fastrelax.get_movemap()
    mp.set_jump(jumps)
    mp.set_bb(bb)
    mp.set_chi(sc)
    fastrelax.set_movemap(mp)
    return fastrelax

def init_rosetta(pose_file, constrain_coords=False, ex1=False, ex2aro=False, ex2=False, really_hard_constraints=False):
    # setup command lines for Rosetta
    opts = [
        "-symmetry:initialize_rigid_body_dofs 1",
        "-use_input_sc",
        "-extrachi_cutoff 1",
        "-unboundrot {}".format(pose_file)
    ]
    if constrain_coords:
        opts.append("-constrain_relax_to_start_coords")  # constrain to starting struture
        opts.append("-relax:coord_constrain_sidechains")  # constrain to starting struture
        opts.append("-relax:ramp_constraints false")  # danger of overpacking structures you have high confidence of
    if really_hard_constraints:
        opts.append("-relax:coord_cst_stdev 0.1")
    if ex1:
        opts.append('-ex1')
    if ex2aro:
        opts.append('-ex2aro')
    if ex2:
        opts.append('-ex2')
    init(extra_options=" ".join(opts))

def symmetric_relax(pose_file, symmetry_file, native_symdef_file=None, output_dir=".", rosetta_out=None, native_file=None, constrain_coords=False, cycles=5, rmsd_map=None):
    # check rmsd_map if parsed
    if rmsd_map is not None:
         rmsd_map = tuple([int(i) if i != "-" else None for i in rmsd_map])

    # assert output directory does exists
    if not Path(output_dir).is_dir():
        raise NotADirectoryError(f"{output_dir} does not exists")

    # set names:
    input_out = str(Path(output_dir).joinpath(Path(pose_file).stem + "_INPUT.pdb"))
    symm_out = str(Path(output_dir).joinpath(Path(pose_file).stem + ".symm"))
    full_out = str(Path(output_dir).joinpath(Path(pose_file).stem + "_full.cif"))
    info_out = str(Path(output_dir).joinpath(Path(pose_file).stem + "_data.csv"))
    if rosetta_out:
        rosetta_out = str(Path(rosetta_out).joinpath(Path(pose_file).stem + ".pdb"))

    # score function and read in symmetric pose
    init_rosetta(pose_file)
    pose = pose_from_file(pose_file)
    pose.conformation().detect_disulfides(vector1_std_pair_unsigned_long_unsigned_long_t())
    if native_file is None:
        native = pose.clone()
    else:
        native = pose_from_file(native_file)
    SetupForSymmetryMover(symmetry_file).apply(pose)
    # fixme: since master merge this does not work with symmetry
    pose.conformation().detect_disulfides(vector1_std_pair_unsigned_long_unsigned_long_t())
    sfxn = ScoreFunctionFactory.create_score_function("ref2015")
    sfxn.score(pose)
    packer = PosePackRotamers(pose, pose_file, "custom")
    pose = packer.pack()
    init_pose = pose.clone()
    if native_symdef_file is None:
        native_symdef_file = symmetry_file
    init_metrics = calculate_metrics(init_pose, native, init_pose, native_symdef_file, symmetry_file, rmsd_map)

    # It can happen that the structure explodes because the output of AF/AFM for
    # instance have bad full structure Rosetta energies.
    # If it is the case we apply constrained relax
    relaxes_done = {"1": False, "2": False, "3": False, "4": False}
    use_old_data = None

    # if constraining to the starting coordinates, only do that and continue
    if constrain_coords:
        fastrelax = make_fastrelax(pose_file, constrain_coords=True, cycles=cycles)
        fastrelax.apply(pose)
        relaxes_done["1"] = True
    # else do a more elaborate relax approach
    else:
        # start out by relaxing the structure according to the input structure if the energy is bad. We dont
        # want to start out just relaxing on the initial structure as it might explode if it has bad energies
        if sfxn.score(pose) > 0:
            fastrelax = make_fastrelax(pose_file, constrain_coords=True, cycles=cycles)
            fastrelax.apply(pose)
            relaxes_done["1"] = True

        # if the energy has improved below 0 we are more confident that the structure will not explode now
        # and we will therefor relax normally
        if sfxn.score(pose) < 0:
            fastrelax = make_fastrelax(pose_file, cycles=cycles)
            fastrelax.apply(pose)
            relaxes_done["2"] = True
        # if the energy has not gone below 0 the pose is turned into the initial pose and then we only optimize sc and jumps
        else:
            pose = init_pose.clone()
            fastrelax = make_fastrelax(pose_file, bb=False, cycles=cycles)
            fastrelax.apply(pose)
            relaxes_done["3"] = True

        # if the interface energy has not improved, then revert back to the initial structure
        if init_metrics["Ienergy"] <= calculate_metrics(pose, native, pose, native_symdef_file, symmetry_file, rmsd_map)["Ienergy"] + 0.1: #0.1 for float precision
            pose = init_pose.clone()
            relaxes_done["4"] = True
            # # FIXME: this is a hack - remove it for END USER!
            # d = pd.read_csv(Path(info_out).parent.parent.parent.joinpath("top_100.csv"))
            # pop, ind, gen = Path(info_out).stem.split(".prepack_")[-1].split("_") # {pop}_{ind}_{gen}
            # use_old_data = d[( (d["pop"] == str(pop)) | (d["pop"] == int(pop)) ) &
            #                  ( (d["ind"] == str(ind)) | (d["ind"] == int(ind))) &
            #                  ( (d["gen"] == str(gen)) | (d["gen"] == int(gen)) ) ]
            # assert len(use_old_data) == 1, "should have 1 unique match!!!!!"

    ####################
    # Now output stuff #
    ####################

    # output the information file
    sfxn.score(pose)
    symmetrysetup = SymmetrySetup()
    symmetrysetup.read_from_file(symmetry_file)
    # name = Path(pose_file).stem + (suffix if suffix else "")
    output_info(pose, native, pose, info_out, native_symdef=native_symdef_file, input_symdef=symmetry_file, relaxes_done=relaxes_done,
                use_old_data= use_old_data, use_rmsd_map=rmsd_map)

    # output the input file and symmetry file
    symmetrysetup.update_dofs_from_pose(pose)
    if CubicSetup(symmetry_file).is_cubic():
        symmetrysetup.make_asymmetric_pose(pose, reset_dofs=True).dump_pdb(input_out)
    else:
        # jumpname from make_symdef_file.pl
        if "JUMP0_to_subunit" in symmetrysetup._dofs.keys():
            symmetrysetup.make_asymmetric_pose(pose, reset_dofs=False).dump_pdb(input_out)
            for i in range(3):
                symmetrysetup._dofs["JUMP0_to_subunit"][i][-1] = None
        else:
            symmetrysetup.make_asymmetric_pose(pose, reset_dofs=True).dump_pdb(input_out)
    symmetrysetup.output(symm_out)

    # output the full symmetrical structure
    pose = pose_from_file(input_out)
    from cubicsym.assembly.cubicassembly import CubicSymmetricAssembly
    cs = CubicSetup(symdef=symm_out)
    cs.make_symmetric_pose(pose)
    CubicSymmetricAssembly.from_pose_input(pose, cs).output(full_out)

    # output the Rosetta structure (IF set)
    if isinstance(rosetta_out, str):
        pose.dump_pdb(rosetta_out)


def main():
    description = ("Wrapper script for the relax protocol in Rosetta. \n" 
                   "In addition to regular relax it fine-tunes the relax protocol to protect it against the input structure blowing up.\n" 
                   "This can happen if the structures produced by AF/AFM are energetically unfavorable when parsed into Rosetta.\n" 
                   "In addition to regular relax outputs it also outputs:\n" 
                   "    1. Monomeric input structure (to be made symmetric in Rosetta with the symmetry file that is also output (see 2.)). Extension is _INPUT.pdb\n"
                   "    2. The symmetry file in which the DOFS (set_dofs lines in the symmetry file) are set to final dofs found in the relax protocol. Extension is .symm\n" 
                   "    3. A fully symmetric structure corresponding to the biological assembly. Extension is _full.cif\n" 
                   "    4. A CSV file containing Iscore/score and Irmsd/rmsd outputs if a --native_file has been parsed. Extension is _data.csv")

    parser = argparse.ArgumentParser(description=description,  formatter_class=RawTextHelpFormatter)
    parser.add_argument("--file", help="Input structure to relax", type=str, required=True)
    parser.add_argument("--symmetry_file", help="Symmetry definition file", type=str, required=True)
    parser.add_argument("--constrain_coords", help="Constrain coordinates only", type=bool, choices=[True, False])
    parser.add_argument("--cycles", help="Relax cycles to perform", type=int, default=5)
    parser.add_argument("--native_file", help="Native file. If not set it will use --file instead (For RMSD calculations).", type=str)
    parser.add_argument("--native_symmetry_file", help="The native symmetry file. If not set it will use --symmetry_file (For RMSD calculations)", type=str)
    parser.add_argument("--output_dir", help="Output directory for all outputs.", type=str, default=".")
    parser.add_argument("--rosetta_out", help="Output the Rosetta symmetric structure.", action="store_true")
    parser.add_argument("--rmsd_map", help="Use an alternative RMSD map than the default", nargs="+")
    args = parser.parse_args()

    symmetric_relax(pose_file=args.file, symmetry_file=args.symmetry_file, native_symdef_file=args.native_symmetry_file,
                    output_dir=args.output_dir,
                    native_file=args.native_file, constrain_coords=args.constrain_coords,
                    cycles=args.cycles, rmsd_map=args.rmsd_map)

if __name__ == "__main__":
    main()

