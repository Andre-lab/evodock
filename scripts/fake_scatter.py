"""New Document"""
#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.pyplot import figure

figure(figsize=(6, 6), dpi=300)


dict_red = {"rosetta": "#bdbdbd", "evodock": "#de2d26", "low": "#3182bd"}


def main():
    x = np.linspace(0, 10, 1)
    y = np.sin(x)
    plt.scatter(x, y, marker="^", color=dict_red["rosetta"], label="RosettaDOCK")
    plt.scatter(
        x,
        y,
        marker="o",
        color=dict_red["evodock"],
        label="EvoDOCK Complete population at final iteration",
    )
    plt.scatter(
        x,
        y,
        marker="o",
        color=dict_red["low"],
        label="EvoDOCK Lowest global energy individuals",
    )

    plt.legend()
    plt.savefig("scatter_with_legend.png")
    plt.close()

    plt.plot(
        x,
        y,
        color=dict_red["rosetta"],
        label="RosettaDOCK",
    )
    plt.plot(
        x,
        y,
        color=dict_red["low"],
        label="EvoDOCK",
    )

    plt.legend(prop={"size": 10})
    plt.savefig("line_with_legend.png")


if __name__ == "__main__":
    main()
