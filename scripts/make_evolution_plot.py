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
    sns.lineplot(x="gen", y="rmsd_from_best", data=df, ax=axes[1])
    fig.tight_layout()
    # fig_name = name.replace(".log", ".png")
    fig_name = "check_evolution.png"
    fig.savefig(fig_name)
    plt.close()
    return fig_name.replace(".png", "")


def create_dataframe(evolution_log_file):
    return pd.read_csv(evolution_log_file)


def main():
    evolution_log_file = sys.argv[-1]

    df = create_dataframe(evolution_log_file)
    plots_during_evolution(df, evolution_log_file)


if __name__ == "__main__":
    main()
