#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ===================== FILES =====================
SN1_NORMAL = "scd_stats_output_normal_sn1.xlsx"
SN1_CANCER = "cancer_scd_stats_output_sn1.xlsx"
SN2_NORMAL = "scd_stats_output_normal_sn2.xlsx"
SN2_CANCER = "cancer_scd_stats_output_sn2.xlsx"

OUT = "Figure_SCD_SN1_SN2_normal_vs_cancer_same_xaxis.pdf"

# ===================== STYLE =====================
FIGSIZE = (7.2, 3.2)

COLOR_NORMAL = "#2e9e59"
COLOR_CANCER = "#e53935"

YLIM = (0.0, 0.27)

# Force both panels to have identical x-axis range and ticks
XLIM = (1.5, 18.5)
XTICKS = np.arange(2, 19, 2)
XTICK_LABELS = [str(int(x)) for x in XTICKS]
# =================================================


def load_scd(path):
    df = pd.read_excel(path)
    return df["Index"], df["mean"], df["std_error"]


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.tick_params(direction="out", length=5, width=1.4)

    # Force same x-axis appearance
    ax.set_xlim(*XLIM)
    ax.set_xticks(XTICKS)
    ax.set_xticklabels(XTICK_LABELS)


def panel_label(ax, label):
    ax.text(
        0.02, 0.95, label,
        transform=ax.transAxes,
        ha="left", va="top",
        fontweight="semibold"
    )


def plot_panel(ax, xN, yN, eN, xC, yC, eC, panel_label_text):
    ax.errorbar(
        xN, yN, yerr=eN,
        fmt="-o",
        markersize=4.0,
        linewidth=1.8,
        capsize=3,
        markerfacecolor=COLOR_NORMAL,
        markeredgecolor=COLOR_NORMAL,
        markeredgewidth=0.5,
        color=COLOR_NORMAL,
        label="Normal"
    )

    ax.errorbar(
        xC, yC, yerr=eC,
        fmt="-^",
        markersize=4.2,
        linewidth=1.8,
        capsize=3,
        markerfacecolor=COLOR_CANCER,
        markeredgecolor=COLOR_CANCER,
        markeredgewidth=0.5,
        color=COLOR_CANCER,
        label="Cancer"
    )

    ax.set_xlabel("Carbon number")
    ax.set_ylim(*YLIM)
    panel_label(ax, panel_label_text)
    style_axes(ax)
    ax.legend(loc="upper right", frameon=False)


def main():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,

        "font.size": 10,
        "axes.labelsize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 1.5,
        "lines.linewidth": 1.8,
        "legend.fontsize": 8,
    })

    x1N, y1N, e1N = load_scd(SN1_NORMAL)
    x1C, y1C, e1C = load_scd(SN1_CANCER)
    x2N, y2N, e2N = load_scd(SN2_NORMAL)
    x2C, y2C, e2C = load_scd(SN2_CANCER)

    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=FIGSIZE,
        sharey=True,
        sharex=True
    )

    plot_panel(ax1, x1N, y1N, e1N, x1C, y1C, e1C, "(a)")
    plot_panel(ax2, x2N, y2N, e2N, x2C, y2C, e2C, "(b)")

    ax1.set_ylabel(r"$S_{CD}$")

    # Apply again after plotting, to make sure Matplotlib does not reset ticks
    for ax in (ax1, ax2):
        ax.set_xlim(*XLIM)
        ax.set_xticks(XTICKS)
        ax.set_xticklabels(XTICK_LABELS)

    fig.tight_layout()
    fig.savefig(OUT, dpi=600, bbox_inches="tight")
    fig.savefig(OUT.replace(".pdf", ".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)

    print("Saved:", Path(OUT).resolve())
    print("Saved:", Path(OUT.replace(".pdf", ".png")).resolve())


if __name__ == "__main__":
    main()
