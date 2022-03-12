#!/usr/bin/env python
# coding: utf-8

import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plots_during_evolution(df, name):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(25, 10))
    ax = sns.lineplot(x="gen", y="best", data=df, ax=axes[0])
    sns.lineplot(x="gen", y="avg", data=df, ax=ax)
    sns.lineplot(x="gen", y="rmsd", data=df, ax=axes[1])
    fig.tight_layout()
    # fig_name = name.replace(".log", ".png")
    fig_name = "check_evolution.png"
    fig.savefig(fig_name)
    plt.close()
    return fig_name.replace(".png", "")


def create_dataframe(evolution_log_file):
    lines = [
        l.strip().split("\t")
        for l in open(evolution_log_file, "r").readlines()
        if "GENERATION" in l
    ]

    data_per_gen = [
        {"gen": int(l[1]), "avg": float(l[2]), "best": float(l[3]), "rmsd": float(l[4])}
        for l in lines
    ]

    return pd.DataFrame(data_per_gen)


def main():
    evolution_log_file = sys.argv[-1]

    df = create_dataframe(evolution_log_file)
    plots_during_evolution(df, evolution_log_file)


if __name__ == "__main__":
    main()
