#!/usr/bin/env python
# coding: utf-8


import glob

import pandas as pd
from pyrosetta import init, pose_from_file
from pyrosetta.rosetta.core.scoring import CA_rmsd, ScoreFunctionFactory

motif_path = "/home/daniel/Workplace/additional_protocol_data-master/motif_dock/xh_16_"


def init_opts():
    opts = [
        # "-mute all",
        "-docking:dock_mcm_first_cycles 1",
        "-docking:dock_mcm_second_cycles 1",
        "-include_current True",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-mh:score:use_ss1 false",
        "-mh:score:use_ss2 false",
        "-mh:score:use_aa1 true",
        "-mh:score:use_aa2 true",
        "-mh:path:scores_BB_BB {}".format(motif_path),
    ]
    return " ".join(opts)


def main():
    pose_input = "./easy_dock/1ZLI.prepack.pdb"
    print(init_opts())
    init(extra_options=init_opts())
    native = pose_from_file(pose_input)
    scfxn_rosetta = ScoreFunctionFactory.create_score_function("ref2015")
    pdbs = [
        pose_from_file(p)
        for p in glob.glob(
            "./ResultsFlexBB_1zli/16092021_Local_1zli_F09CR09_0/1zli/*/*pdb"
        )
    ]

    results = []
    for p in pdbs:
        score = scfxn_rosetta.score(p)
        rmsd = CA_rmsd(native, p)
        print("{} {}".format(score, rmsd))
        results.append({"score": score, "rmsd": rmsd})

    df = pd.DataFrame(results)
    ax = df.plot.scatter(x="rmsd", y="score")
    ax.get_figure().savefig("1zli_bound_local_docking.png")


if __name__ == "__main__":
    main()
