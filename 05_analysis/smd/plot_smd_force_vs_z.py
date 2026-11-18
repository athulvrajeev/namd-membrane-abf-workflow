#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================

DATA_FILE = Path("bare_cancer_FvsZ_240ns.dat")

OUT = "Figure_SMD_bare_cancer_FvsZ_240.pdf"

TITLE = "Bare AuNP — cancer-associated membrane model"

XLABEL = r"Reaction coordinate, $z$ ($\mathrm{\AA}$)"
YLABEL = "SMD force"

# If your force is already in pN, use:
# YLABEL = "SMD force (pN)"

FIGSIZE = (7.2, 3.4)

# Binning/smoothing settings
BIN_WIDTH = 0.5      # Å
ROLLING_WINDOW = 5   # number of bins for smoothing

PLOT_RAW = True

# ============================================================


def set_pccp_style():
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",

        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,

        "font.size": 10,
        "axes.labelsize": 12,
        "axes.titlesize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 1.5,

        "legend.fontsize": 8,
    })


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.tick_params(direction="out", length=5, width=1.4)
    ax.axhline(0, color="0.5", linestyle=":", linewidth=1.0, alpha=0.8)
    ax.axvline(0, color="0.5", linestyle="--", linewidth=1.0, alpha=0.8)


def read_smd_file(fname: Path):
    data = np.loadtxt(fname)

    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("Input file must contain at least two columns: z and force.")

    z = data[:, 0]
    force = data[:, 1]

    mask = np.isfinite(z) & np.isfinite(force)
    return z[mask], force[mask]


def bin_average(x, y, bin_width):
    xmin, xmax = np.min(x), np.max(x)
    bins = np.arange(xmin, xmax + bin_width, bin_width)

    centers = 0.5 * (bins[:-1] + bins[1:])
    means = np.full_like(centers, np.nan, dtype=float)
    stds = np.full_like(centers, np.nan, dtype=float)

    for i in range(len(centers)):
        mask = (x >= bins[i]) & (x < bins[i + 1])

        if np.any(mask):
            means[i] = np.mean(y[mask])
            stds[i] = np.std(y[mask])

    valid = np.isfinite(means)
    return centers[valid], means[valid], stds[valid]


def rolling_average(y, window):
    if window <= 1:
        return y

    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")


def main():
    set_pccp_style()

    z, force = read_smd_file(DATA_FILE)

    # Sort by reaction coordinate for clean Force-vs-Z plotting
    order = np.argsort(z)
    z_sorted = z[order]
    force_sorted = force[order]

    z_bin, force_mean, force_sd = bin_average(z_sorted, force_sorted, BIN_WIDTH)
    force_smooth = rolling_average(force_mean, ROLLING_WINDOW)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    if PLOT_RAW:
        ax.plot(
            z_sorted,
            force_sorted,
            color="0.75",
            linewidth=0.6,
            alpha=0.45,
            label="Raw data"
        )

    ax.plot(
        z_bin,
        force_smooth,
        color="#e53935",
        linewidth=2.2,
        label="Smoothed profile"
    )

    ax.fill_between(
        z_bin,
        force_mean - force_sd,
        force_mean + force_sd,
        color="#e53935",
        alpha=0.15,
        linewidth=0,
        label="Local SD"
    )

    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_title(TITLE)

    style_axes(ax)

    ax.legend(frameon=False, loc="best")

    fig.tight_layout()

    fig.savefig(OUT)
    fig.savefig(OUT.replace(".pdf", ".png"), dpi=600)

    print(f"Saved: {OUT}")
    print(f"Saved: {OUT.replace('.pdf', '.png')}")


if __name__ == "__main__":
    main()
