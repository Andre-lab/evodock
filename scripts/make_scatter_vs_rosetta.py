#!/usr/bin/env python
# coding: utf-8

import glob
import os
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def logfile_to_dataframe(log):
    idfile = log.replace(".log", "")
    with open(log, "r") as f:
        lines = [l.strip() for l in f.readlines() if len(l) > 0]

    final_lines = lines[-2:]
    scores = np.array(final_lines[0].split(","), dtype=float)
    rmsds = np.array(final_lines[1].split(","), dtype=float)
    individuals = [
        {"id": idfile, "score": x[0], "rmsd": x[1]} for x in list(zip(scores, rmsds))
    ]

    df = pd.DataFrame(individuals)
    df["source"] = ["evodock"] * len(df)
    # print(df)
    return df


def best_pdb(df):
    df = df[df.source == "evodock"]
    min_value = df[df.rmsd == df.rmsd.min()]
    print(min_value.head())
    print(min_value.iloc[0].id)
    pdb = (
        "/".join(min_value.iloc[0].id.split("/")[:-1])
        + "/evolution_f03_cr09_final_docked_evo.pdb"
    )
    if os.path.isfile(pdb):
        shutil.copy(pdb, "best_found.pdb")


def print_plot(name, df, interface=False, zoom=False):
    ax = sns.scatterplot(
        x="rmsd", y="score", data=df, hue="source", alpha=0.4, markers=["s"],
    )
    if zoom:
        ax.set_xlim([0, min(df["rmsd"].values.tolist()) + 10])
        ax.set_ylim(
            [
                min(df["score"].values.tolist()) - 10,
                min(df["score"].values.tolist()) + 30,
            ]
        )
    if interface:
        ax.set_xlabel("iRMSD")
        ax.set_ylabel("interface REU (Rosetta Energy Units)")
    else:
        ax.set_xlabel("RMSD")
        ax.set_ylabel("REU (Rosetta Energy Units)")
    fig = ax.get_figure()
    fig.tight_layout()

    name = "interface_" + name if interface else name
    name = "zoom_" + name if zoom else name
    fig.savefig(name)
    plt.close()


def main():
    input_files = glob.glob(sys.argv[-1])
    print(input_files)
    results = []
    for ifile in input_files:
        results.append(logfile_to_dataframe(ifile))

    df_rosetta = pd.read_csv(
        "~/projects/RosettaDockingData/LOCAL_DATASETS/rosetta_1ZLI.csv"
    )
    print(df_rosetta.head())
    if "interface" in input_files[0]:
        df_rosetta = df_rosetta[["I_sc", "Irms", "description"]]
        df_rosetta.columns = ["score", "rmsd", "id"]
    else:
        df_rosetta = df_rosetta[["score", "rms", "description"]]
        df_rosetta.columns = ["score", "rmsd", "id"]

    df_rosetta["source"] = ["rosetta"] * len(df_rosetta)
    df = pd.concat(results, ignore_index=True)
    print(df.head())
    df = pd.concat([df_rosetta, df], ignore_index=True)
    print(df_rosetta.head())

    name = input_files[0].split("/")[0] + ".png"
    interface = "interface" in input_files[0]
    print_plot(name, df, interface)
    print_plot(name, df, interface, zoom=True)
    best_pdb(df)


if __name__ == "__main__":
    main()
