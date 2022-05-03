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

dict_red = {"rosetta": "#bdbdbd", "evodock": "#de2d26", "low": "#3182bd"}


def colored_scatter(x, y, c=None, m=None):
    def scatter(*args, **kwargs):
        args = (x, y)
        if c is not None:
            kwargs["c"] = c
        if c == "orange":
            scatter_alpha = 0.1
        else:
            scatter_alpha = 0.5

        kwargs["alpha"] = scatter_alpha
        kwargs["marker"] = m
        plt.scatter(*args, **kwargs)

    return scatter


def logfile_to_dataframe(log, interface=False):
    idfile = log.replace(".log", "")
    with open(log, "r") as f:
        lines = [l.strip() for l in f.readlines() if len(l) > 0]

    d = log.split("/")
    d[-1] = d[-1].replace("popul_", "interface_")
    print("/".join(d))
    with open("/".join(d), "r") as f:
        interface_lines = [l.strip() for l in f.readlines() if len(l) > 0]

    if len(lines) > 1:
        final_lines = lines[-2:]
        final_interface_lines = interface_lines[-2:]
        scores = np.array(final_lines[0].split(","), dtype=float)
        rmsds = np.array(final_lines[1].split(","), dtype=float)

        score_interface = np.array(final_interface_lines[0].split(","),
                                   dtype=float)
        rmsds_interface = np.array(final_interface_lines[1].split(","),
                                   dtype=float)

        # individuals = [
        #     {"id": idfile, "score": x[0], "rmsd": x[1]} for x in list(zip(scores, rmsds))
        # ]

        r = range(min(len(scores), len(score_interface)))
        print(rmsds_interface[:10])
        print(rmsds[:10])

        individuals = [{
            "id": idfile,
            "score": scores[x],
            "rmsd": rmsds[x],
            "I_sc": score_interface[x],
            "Irms": rmsds_interface[x]
        } for x in r]

        df = pd.DataFrame(individuals)
        df["source"] = ["evodock"] * len(df)
        df_low = df[df.score == df.score.min()]
        df_low.source = "lowest global energy individuals"
        df = pd.concat([df_low, df])
        if interface:
            df["rmsd"] = df["Irms"]
            # df["score"] = df["I_sc"]

        return df
    else:
        return pd.DataFrame([])


def best_pdb(df):
    df = df[df.source == "evodock"]
    min_value = df[df.rmsd == df.rmsd.min()]
    print(min_value.head())
    print(min_value.iloc[0].id)
    pdb = ("/".join(min_value.iloc[0].id.split("/")[:-1]) +
           "/evolution_f03_cr09_final_docked_evo.pdb")
    if os.path.isfile(pdb):
        shutil.copy(pdb, "best_found.pdb")


def print_plot(name, df, interface=False, zoom=False):
    print("evodock energies: ({}, {})".format(
        df[df.source == "evodock"].score.min(),
        df[df.source == "evodock"].score.max()))
    print("rosetta energies: ({}, {})".format(
        df[df.source == "rosetta"].score.min(),
        df[df.source == "rosetta"].score.max()))

    print("evodock Irmsd: ({}, {})".format(
        df[df.source == "evodock"].Irms.min(),
        df[df.source == "evodock"].Irms.max()))
    print("rosetta Irmsd: ({}, {})".format(
        df[df.source == "rosetta"].rmsd.min(),
        df[df.source == "rosetta"].rmsd.max()))

    min_evodock = df[df.source == "evodock"].score.min()
    print("lowest energy evodock {} ({})".format(
        df[df.score == min_evodock].iloc[0].score,
        df[df.score == min_evodock].iloc[0].Irms))
    min_rosetta = df[df.source == "rosetta"].score.min()
    print("lowest energy rosetta {} ({})".format(
        df[df.score == min_rosetta].iloc[0].score,
        df[df.score == min_rosetta].iloc[0].rmsd))

    ax = sns.scatterplot(
        x="rmsd",
        y="score",
        data=df,
        hue="source",
        alpha=0.4,
        markers=["s"],
    )

    colors = {
        "rosetta": "#bdbdbd",
        "evodock": "#de2d26",
        "lowest global energy individuals": "#3182bd"
    }
    markers = {
        "rosetta": "^",
        "evodock": "o",
        "lowest global energy individuals": "o"
    }
    names_sort = ["rosetta", "evodock", "lowest global energy individuals"]
    legends = []

    col_x, col_y = "rmsd", "score"

    g = sns.JointGrid(x=col_x, y=col_y, data=df)

    for name in names_sort:
        df_group = df[df.source == name]
        legends.append(name)
        color = colors[name]
        marker = markers[name]
        g.plot_joint(
            colored_scatter(df_group[col_x], df_group[col_y], color, marker), )
        sns.distplot(
            df_group[col_x].values,
            ax=g.ax_marg_x,
            color=color,
        )
        sns.distplot(df_group[col_y].values,
                     ax=g.ax_marg_y,
                     color=color,
                     vertical=True)

    ax = g.ax_joint

    # lims = (2, 30) if interface else (10, 30)
    lims = (2, 60) if interface else (10, 30)
    if zoom:
        ax.set_xlim([0, min(df["rmsd"].values.tolist()) + 10])
        ax.set_ylim([
            min(df["score"].values.tolist()) - lims[0],
            min(df["score"].values.tolist()) + lims[1],
        ])
    props = dict(boxstyle="round", facecolor="none", alpha=0.0)

    # plt.text(
    #     1.10,
    #     1.10,
    #     "1zli",
    #     horizontalalignment="center",
    #     weight="bold",
    #     va="center",
    #     fontsize=18,
    #     transform=ax.transAxes,
    #     bbox=props,
    # )

    if interface:
        ax.set_xlabel("iRMSD")
        #ax.set_ylabel("interface REU (Rosetta Energy Units)")
        ax.set_ylabel("REU (Rosetta Energy Units)")
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
    interface = True
    for ifile in input_files:
        results.append(logfile_to_dataframe(ifile, interface))

    df_rosetta = pd.read_csv(
        "~/projects/RosettaDockingData/LOCAL_DATASETS/rosetta_1ZLI.csv")
    print(df_rosetta.head())
    if interface:
        # df_rosetta = df_rosetta[["I_sc", "Irms", "description"]]
        df_rosetta = df_rosetta[["score", "Irms", "description"]]
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
    # interface = "interface" in input_files[0]
    interface = True
    print_plot(name, df, interface)
    print_plot(name, df, interface, zoom=True)
    # best_pdb(df)


if __name__ == "__main__":
    main()
