#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ===================== USER SETTINGS =====================
FIGSIZE = (7.2, 3.2)
OUT = "Figure_density_normal_vs_cancer_final.pdf"

SMOOTH = True

SAVGOL_WINDOW = 7
SAVGOL_POLY = 2

TAIL_EPS = 5e-5

normal_files = {
    "POPC": "POPC_density_profile.csv",
    "POPE": "POPE_density_profile.csv",
    "POPS": "POPS_density_profile.csv",
    "PSM":  "PSM_density_profile.csv",
    "CHL":  "CHL1_density_profile.csv",
}

cancer_files = {
    "POPC": "can_POPC_density_profile.csv",
    "POPE": "can_POPE_density_profile.csv",
    "POPS": "can_POPS_density_profile.csv",
    "PSM":  "can_PSM_density_profile.csv",
    "CHL":  "can_CHL1_density_profile.csv",
}

colors = {
    "POPC": "#4d4d4d",
    "POPE": "#e53935",
    "POPS": "#1565c0",
    "PSM":  "#2e9e59",
    "CHL":  "#a970d6",
}

markers = {
    "POPC": "o",
    "POPE": "^",
    "POPS": "s",
    "PSM":  "D",
    "CHL":  "p",
}

MARK_EVERY = 5

XLIM = (-30, 30)
YLIM = (0.0, 0.035)
# =========================================================


def read_density_csv(path: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path)
    x = df["Z-coordinate"].to_numpy()
    y = df["Density"].to_numpy()
    return x, y


def safe_smooth(y: np.ndarray, window: int, poly: int) -> np.ndarray:
    if window % 2 == 0:
        window += 1

    window = max(window, poly + 3)

    try:
        from scipy.signal import savgol_filter
        return savgol_filter(
            y,
            window_length=window,
            polyorder=poly,
            mode="interp"
        )
    except Exception:
        kernel = np.ones(window) / window
        return np.convolve(y, kernel, mode="same")


def remove_artifacts(
    y_raw: np.ndarray,
    y_proc: np.ndarray,
    tail_eps: float
) -> np.ndarray:
    y = np.clip(y_proc, 0.0, None)
    y[y_raw < tail_eps] = 0.0
    return y


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.tick_params(direction="out", length=5, width=1.4)


def panel_label(ax, label: str):
    ax.text(
        0.02,
        0.95,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontweight="semibold",
    )


def main():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,

        "font.size": 10,
        "axes.labelsize": 12,
        "axes.titlesize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 1.5,
        "lines.linewidth": 1.8,
        "legend.fontsize": 8,
    })

    fig, axes = plt.subplots(
        1,
        2,
        figsize=FIGSIZE,
        sharex=True,
        sharey=True
    )

    ax1, ax2 = axes

    # ---------- Panel (a): Normal ----------
    for comp, path in normal_files.items():
        x, y_raw = read_density_csv(path)

        if SMOOTH:
            y = safe_smooth(y_raw, SAVGOL_WINDOW, SAVGOL_POLY)
        else:
            y = y_raw.copy()

        y = remove_artifacts(y_raw, y, TAIL_EPS)

        ax1.plot(
            x, y,
            color=colors[comp],
            linestyle="-",
            linewidth=1.8,
            marker=markers[comp],
            markersize=3.5,
            markerfacecolor=colors[comp],
            markeredgecolor=colors[comp],
            markeredgewidth=0.4,
            markevery=MARK_EVERY,
            label=comp
        )

    ax1.set_xlabel("Distance from membrane centre (Å)")
    ax1.set_ylabel(r"Number density (atoms Å$^{-3}$)")
    ax1.set_xlim(*XLIM)
    ax1.set_ylim(*YLIM)
    panel_label(ax1, "(a)")
    style_axes(ax1)

    # ---------- Panel (b): Cancer ----------
    for comp, path in cancer_files.items():
        x, y_raw = read_density_csv(path)

        if SMOOTH:
            y = safe_smooth(y_raw, SAVGOL_WINDOW, SAVGOL_POLY)
        else:
            y = y_raw.copy()

        y = remove_artifacts(y_raw, y, TAIL_EPS)

        ax2.plot(
            x, y,
            color=colors[comp],
            linestyle="-",
            linewidth=1.8,
            marker=markers[comp],
            markersize=3.5,
            markerfacecolor=colors[comp],
            markeredgecolor=colors[comp],
            markeredgewidth=0.4,
            markevery=MARK_EVERY
        )

    ax2.set_xlabel("Distance from membrane centre (Å)")
    panel_label(ax2, "(b)")
    style_axes(ax2)

    handles, labels = ax1.get_legend_handles_labels()

    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02)
    )

    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(OUT, dpi=600)
    plt.close(fig)

    print(f"Saved: {Path(OUT).resolve()}")


if __name__ == "__main__":
    main()
