#!/usr/bin/env python
# coding: utf-8

import glob
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
    # print(df)
    return df


def main():
    input_files = glob.glob(sys.argv[-1])

    print(input_files)
    results = []
    for ifile in input_files:
        results.append(logfile_to_dataframe(ifile))

    df = pd.concat(results, ignore_index=True)
    ax = sns.scatterplot(
        x="rmsd", y="score", data=df, hue="id", alpha=0.4, markers=["s"],
    )
    ax.set_xlim([0, max(df["rmsd"].values.tolist())])
    ax.get_legend().remove()
    if "interface" in input_files[0]:
        ax.set_xlabel("iRMSD")
        ax.set_ylabel("interface REU (Rosetta Energy Units)")
    else:
        ax.set_xlabel("RMSD")
        ax.set_ylabel("REU (Rosetta Energy Units)")
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig("scatterplot.png")
    plt.close()

    print("min rmsd {}".format(df.rmsd.min()))


main()
