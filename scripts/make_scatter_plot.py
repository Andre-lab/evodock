#!/usr/bin/env python
# coding: utf-8

import glob
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def logfile_to_dataframe(log):
    idfile = log.replace(".csv", "")
    data = pd.read_csv(log).iloc[-1]
    score_columns = [c for c in data.index if c.startswith("sc_")]
    rmsd_columns = [c for c in data.index if c.startswith("rmsd_")]
    isc_columns = [c for c in data.index if c.startswith("Isc_")]
    irmsd_columns = [c for c in data.index if c.startswith("Irmsd_")]
    scores = data[score_columns].tolist()
    rmsds = data[rmsd_columns].tolist()
    iscs = data[isc_columns].tolist()
    irmsds = data[irmsd_columns].tolist()
    individuals = {
        "id": [idfile] * len(scores),
        "score": scores,
        "rmsd": rmsds,
        "Isc": iscs,
        "Irmsd": irmsds,
    }

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
        x="rmsd",
        y="score",
        data=df,
        hue="id",
        alpha=0.4,
        markers=["s"],
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
