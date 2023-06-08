#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run relax with symmetry and allows one to output the final symmetry file and the input file centered
@Author: Mads Jeppesen
@Date: 6/15/22
"""
import argparse
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta import init, pose_from_file, standard_task_factory, Pose
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, IncludeCurrent, InitializeFromCommandline
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from symmetryhandler.symmetrysetup import SymmetrySetup
from pathlib import Path
import pandas as pd
from src.dock_metric import CubicDockMetric
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_jump_num
from cubicsym.cubicsetup import CubicSetup
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.pose.symmetry import extract_asymmetric_unit

# fixme it does not like this relative import so will just add this here (copy from it!)
# from ..scripts.prepacking import PosePackRotamers
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent,
    InitializeFromCommandline,
    NoRepackDisulfides,
    RestrictToRepacking,
)
from pyrosetta import get_fa_scorefxn
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

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

    # def run(self):
    #     pack_pdb = self.pack()
    #     pack_pdb.dump_pdb(self.output_name)
    #     return pack_pdb
    #
    # def pack(self):
    #     if self.mode == "standard":
    #         return PoseStandardPacker(self.pose, self.output_name).pack()
    #     if self.mode == "custom":
    #         return PosePackRotamers(self.pose, self.output_name).pack()

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


import numpy as np

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


# def get_score_terms(pose):
#     # buf = ostringstream()
#     # # THESE ARE UNWEIGHTHED:
#     # # from https://graylab.jhu.edu/pyrosetta/downloads/documentation/pyrosetta4_online_format/PyRosetta4_Workshop3_Scoring.pdf
#     # # Unweighted, individual component energies of each residue in a structure are stored in the Pose object and can be accessed by its energies() object.
#     # # For example, to break the energy down into each residue’s contribution use:
#     # # print ras.energies().show(<n>)
#     # pose.energies().show(buf)
#     #return {i.split()[0]: i.split()[1] for i in [l.replace("total_energy", "") for l in buf.str().split("\n") if "total_energy" in l]}
#     return {(pose.scores.keys())}

def calculate_metrics(pose, native, input_pose, native_symdef, input_symdef):
    # get the translational dfo
    jumpints, dofints, trans_mags = [], [], []
    foldid =CubicSetup.get_jumpidentifier_from_pose(pose)
    for jump, dof, trans in zip([f"JUMP{foldid}fold1", f"JUMP{foldid}fold111"], [ "z", "x"], [2000, 1000]):
        jumpints.append(sym_dof_jump_num(pose, jump))
        dofints.append(dof_map[dof])
        trans_mags.append(trans)
    info = dict(pose.scores.items()) # all score terms
    dockmetric = CubicDockMetric(native, input_pose, native_symdef=native_symdef, input_symdef=input_symdef,
                                 jump_ids=jumpints, dof_ids=dofints, trans_mags=trans_mags)
    info["energy"] = info.pop("total_score")
    info["rmsd"] = dockmetric.ca_rmsd(pose)
    info["Irmsd"] = dockmetric.interface_rmsd(pose)
    info["Ienergy"] = dockmetric.interaction_energy(pose)
    return info

def outpout_info(pose, native, input_pose ,outpath, native_symdef, input_symdef, relaxes_done, use_old_data = None):
    # get the translational dfos
    info = calculate_metrics(pose, native, input_pose, native_symdef, input_symdef)
    if use_old_data is not None:
        # owverwrite
        info = {k: np.NaN for k, v in info.items()}
        info["rmsd"] = use_old_data["rmsd"].values[0]
        info["Irmsd"] = use_old_data["Irmsd"].values[0]
        info["energy"] = use_old_data["score"].values[0]
        info["Ienergy"] = use_old_data["Iscore"].values[0]
    for k, v in relaxes_done.items():
        info[k] = v
    if Path(outpath).is_dir():
        outpath = f"{outpath}/info.csv"
    pd.DataFrame(info, index=[0]).to_csv(Path(outpath), index=False)

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

def symmetric_relax(pose_file, symmetry_file, native_symdef_file, cycles=5, constrain_coords=False,
                    rosetta_out=".", input_out=".", symm_out=".", ex1=False, ex2aro=False, ex2=False, info_out=None, native_file=None, suffix=None,
                    bb=True, sc=True, jumps=False):

    # score function and read in symmetric pose
    init_rosetta(pose_file)
    pose = pose_from_file(pose_file)
    pose.conformation().detect_disulfides()
    SetupForSymmetryMover(symmetry_file).apply(pose)
    native = pose_from_file(native_file)
    #native.conformation().detect_disulfides()
    pose.conformation().detect_disulfides()
    sfxn = ScoreFunctionFactory.create_score_function("ref2015")
    sfxn.score(pose)
    # init_pose.conformation().detect_disulfides()
    #
    packer = PosePackRotamers(pose, pose_file, "custom")
    pose = packer.pack()
    #
    init_pose = pose.clone()
    init_metrics = calculate_metrics(init_pose, native, init_pose, native_symdef_file, symmetry_file)


    # todo: if this does not work.
    # could do expensive and might not work, do interface refinement.
    # Check the interface energy, and if not below the current value,
    # what didnt work:
        # extensive constrained relax (15 cycles)
        # constrained relax followed by

    # It can happen that the structure explodes because the output of AF/AFM for
    # instance have bad full structure Rosetta energies.
    # If it is the case we apply constrained relax
    relaxes_done = {"1": False, "2": False, "3": False, "4": False}

    if sfxn.score(pose) > 0:
        fastrelax = make_fastrelax(pose_file, constrain_coords=True)
        fastrelax.apply(pose)
        relaxes_done["1"] = True
    # if the energy has improved below 0 we are more confident that the structure will not explode now
    # and we will therefor relax normally
    if sfxn.score(pose) < 0:
        fastrelax = make_fastrelax(pose_file)
        fastrelax.apply(pose)
        relaxes_done["2"] = True
    # if the energy has not gone below 0 the pose is turned into the initial pose and then we only optimize sc and jumps
    else:
        pose = init_pose.clone()
        fastrelax = make_fastrelax(pose_file, bb=False)
        fastrelax.apply(pose)
        relaxes_done["3"] = True

    # if the interface energy has not improved, then revert back to the initial structure
    use_old_data = None
    if init_metrics["Ienergy"] <= calculate_metrics(pose, native, pose, native_symdef_file, symmetry_file)["Ienergy"] + 0.1: #0.1 for float precision
        pose = init_pose.clone()
        relaxes_done["4"] = True
        # FIXME: this is a hack - remove it for END USER!
        d = pd.read_csv(Path(info_out).parent.parent.parent.joinpath("top_100.csv"))
        pop, ind, gen = Path(info_out).stem.split(".prepack_")[-1].split("_") # {pop}_{ind}_{gen}
        use_old_data = d[( (d["pop"] == str(pop)) | (d["pop"] == int(pop)) ) &
                         ( (d["ind"] == str(ind)) | (d["ind"] == int(ind))) &
                         ( (d["gen"] == str(gen)) | (d["gen"] == int(gen)) ) ]
        assert len(use_old_data) == 1, "should have 1 unique match!!!!!"

    # Sometimes structures can explode from relax. Therefore we run 1 cycle and check the RMSD before moving on
    # pose_init = pose.clone()
    # init_score = sfxn.score(pose)
    # fastrelax_init = make_fastrelax(1, jumps, bb, sc)
    # fastrelax_init.apply(pose_init)
    # pose_init_asym = Pose()
    # pose_asym = Pose()
    # extract_asymmetric_unit(pose, pose_asym, False)
    # extract_asymmetric_unit(pose_init, pose_init_asym, False)
    # if CA_rmsd(pose_init, pose_asym) > 10:
    #     after_score = sfxn.score(pose_init)
    #     print(f"The structure is about to explode!")
    #     print(f"Energy was initally {init_score} and now the current RMSD to the input structure is {CA_rmsd(pose_init, pose_asym)}")
    #     print(f"A single cycle with constrained relax will be carried out without moving the interface")
    #     init_rosetta(True, ex1, ex2aro, ex2)
    #     pose_init = pose.clone()
    #     fastrelax_constrain = make_fastrelax(1, False, bb, sc)
    #     fastrelax_constrain.apply(pose_init)
    #     pose_asym = Pose()
    #     extract_asymmetric_unit(pose.clone(), pose_asym, False)
    #     # try again with really hard constraints with 2 cycles if the RMSD is more than 5 Å!
    #     if CA_rmsd(pose_init, pose_asym) > 5:
    #         init_rosetta(True, ex1, ex2aro, ex2, really_hard_constraints=True)
    #         pose_init = pose.clone()
    #         fastrelax_constrain = make_fastrelax(2, False, bb, sc)
    #         fastrelax_constrain.apply(pose)
    #         assert CA_rmsd(pose_init,pose_asym) < 5, "Constrained relax could not save it. Relax does not give meaningfull results on this structure!"
    #     after_score = sfxn.score(pose_init)
    #     pose.assign(pose_init)
    #     print(f"With constrained relax the energy is now: {after_score} with RMSD to the input structure {CA_rmsd(pose_init, pose_asym)}")
    #     init_rosetta(constrain_coords, ex1, ex2aro, ex2)


    # Now output the final symmetry file
    sfxn.score(pose)
    symmetrysetup = SymmetrySetup()
    symmetrysetup.read_from_file(symmetry_file)
    name = Path(pose_file).stem + (suffix if suffix else "")
    outpout_info(pose, native, pose, info_out, native_symdef=native_symdef_file, input_symdef=symmetry_file, relaxes_done=relaxes_done,
                 use_old_data= use_old_data)
    pose.dump_pdb(rosetta_out + f"/rosetta_{name}.pdb") # full out
    symmetrysetup.update_dofs_from_pose(pose)
    # FIXME: dont_rest is a HACK FOR CUBIC FOR NOW
    # symmetrysetup.make_asymmetric_pose(pose, reset_dofs=True, dont_reset=["JUMPHFfold1111_subunit"]).dump_pdb(input_out + f"/input_{name}.pdb") # input out
    symmetrysetup.make_asymmetric_pose(pose, reset_dofs=True).dump_pdb(input_out + f"/input_{name}.pdb") # input out
    symmetrysetup.output(symm_out + f"/{name}.symm") # sym out
    # make a symmetric pose again with the outputs
    # open the symmetryfile again and remove this line:  JUMPHFfold1111_subunit
    # with open(symm_out + f"/{name}.symm", "r") as f:
    #     lines = f.readlines()
    # with open(symm_out + f"/{name}.symm", "w") as f:
    #     for line in lines:
    #         if not "set_dof JUMPHFfold1111_subunit" in line:
    #             f.write(line)
    pose = pose_from_file(input_out + f"/input_{name}.pdb")
    SetupForSymmetryMover(symm_out + f"/{name}.symm").apply(pose)
    pose.dump_pdb(rosetta_out + f"/rosetta_recreated_from_input_and_symm_file_{name}.pdb")

def main():
    description = "Wrapper script for the relax binary. In addition to regular relax which outputs the full symmetric structure, this scripts " \
                  "also outputs an input structure (to be made symmetric in Rosetta) and the symmetry file which dofs values (in set_dofs) " \
                  "are set to the end point of the relax protocol and a csv file of the final unweighted score terms, total score and RMSD." \
                  "a recreated file using the outputtet input file and symmetry file is also created."

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--file", help="pose to relax", type=str, required=False)
    parser.add_argument("--symmetry_file", help="symmetry definition file", type=str)
    parser.add_argument("--native_file", help="native file. If not set it will use --file instead.", type=str)
    parser.add_argument("--native_symmetry_file", help="The native symmetry file", type=str)
    parser.add_argument("--constrain_coords", help="Constrain to starting coordinates", action="store_true")
    parser.add_argument("--rosetta_out", help="output of the full symmetric structure", type=str, default=".")
    parser.add_argument("--input_out", help="output of the input structure", type=str, default=".")
    parser.add_argument("--sym_out", help="output of the symmetry file", type=str, default=".")
    parser.add_argument("--info_out", help="output of the csv file containing score terms and RMSD", type=str, default=".")
    parser.add_argument("--suffix", help="suffix given to all output structures", type=str)
    # fixme: need to do packing before and after, and to be even more thorough prob do relax as well
    # parser.add_argument("--move_jumpdofs", help="When calculating the interaction energy move these jumpdofs with a certain amount. "
    #                                             "A jumpdof and is pertubation is specified by its jump name, dof name and an amount as follows: "
    #                                             "<jumpname>:<dofname>:<amount>. Example: Jump1:z:1000.", nargs="+", required=True)
    # minization options
    parser.add_argument("--bb", help="minimize backbone.", action="store_true")
    parser.add_argument("--sc", help="minimize sidecahins.", action="store_true")
    parser.add_argument("--jumps", help="minimize jumps.", action="store_true")
    # direct relax optionsoptions
    parser.add_argument("--cycles", help="cycles to use.", type=int, default=15)
    parser.add_argument("--ex1", help="Do ex1 rotamer sampling", action="store_true")
    parser.add_argument("--ex2aro", help="Do ex2aro rotamer sampling", action="store_true")
    parser.add_argument("--ex2", help="Do ex2 rotamer sampling", action="store_true")
    args = parser.parse_args()

    symmetric_relax(args.file, args.symmetry_file, args.native_symmetry_file, args.cycles, args.constrain_coords,
                    args.rosetta_out, args.input_out, args.sym_out, args.ex1, args.ex2aro, args.ex2, args.info_out, args.native_file, args.suffix,
                     args.bb, args.sc, args.jumps)

if __name__ == "__main__":

    main()

