#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ===================== USER SETTINGS =====================

# Use your rightmost-normalized PMF files here
PMF_NORMAL = "merge_nor_gold_right0_normalized.pmf"
PMF_CANCER = "merge_can_gold_right0_normalized.pmf"

OUT = "Figure_PMF_bare_AuNP_normal_vs_cancer_no_smoothing.pdf"

# PCCP-style two-column figure size
FIGSIZE = (7.2, 3.2)

COLOR_NORMAL = "#2e9e59"  # green
COLOR_CANCER = "#e53935"  # red

XLABEL = "Distance from membrane COM (Å)"
YLABEL = "PMF (kcal mol$^{-1}$)"

# Optional axis limits
XLIM = (-35, 35)

# Use None for automatic y-axis limits
YLIM_NORMAL = None
YLIM_CANCER = None

# Save the plotted/referenced data
SAVE_REFERENCED_FILES = True
OUT_REF_NORMAL = "merge_nor_gold_right0_plotted_raw.dat"
OUT_REF_CANCER = "merge_can_gold_right0_plotted_raw.dat"

# Safety option:
# Even if the input PMF is already normalized, this ensures the plotted
# rightmost point is exactly 0 kcal/mol.
REANCHOR_RIGHTMOST = True

# =========================================================


def load_xy(path: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Read first two numeric columns from a PMF file.
    Column 1 = reaction coordinate
    Column 2 = PMF
    Comment/header lines beginning with #, @, or ; are ignored.
    """
    rows = []

    with open(path, "r") as f:
        for line in f:
            s = line.strip()

            if not s:
                continue

            if s.startswith(("#", "@", ";")):
                continue

            parts = s.replace(",", " ").split()

            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    pass

            if len(vals) >= 2:
                rows.append((vals[0], vals[1]))

    if not rows:
        raise ValueError(f"No numeric PMF data found in: {path}")

    arr = np.array(rows, dtype=float)
    x = arr[:, 0]
    y = arr[:, 1]

    order = np.argsort(x)
    return x[order], y[order]


def anchor_to_rightmost_x(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, float]:
    """
    Shift PMF so that the rightmost x-value has PMF = 0.
    This only changes the arbitrary additive constant of the PMF.
    """
    idx = np.argmax(x)
    ref_x = float(x[idx])
    offset = float(y[idx])

    y_ref = y - offset

    return y_ref, offset, ref_x


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
        fontweight="bold",
    )


def main():
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
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 1.5,

        "lines.linewidth": 2.6,
        "legend.fontsize": 10,
    })

    # ---- Load normalized PMF files ----
    xN, yN = load_xy(PMF_NORMAL)
    xC, yC = load_xy(PMF_CANCER)

    # ---- Optional safety re-anchoring at rightmost x-value ----
    if REANCHOR_RIGHTMOST:
        yN, offN, refN = anchor_to_rightmost_x(xN, yN)
        yC, offC, refC = anchor_to_rightmost_x(xC, yC)
    else:
        offN, refN = 0.0, float(xN[np.argmax(xN)])
        offC, refC = 0.0, float(xC[np.argmax(xC)])

    # ---- Save plotted data ----
    if SAVE_REFERENCED_FILES:
        np.savetxt(
            OUT_REF_NORMAL,
            np.column_stack([xN, yN]),
            header=(
                "x(Å)  PMF_rightmost_referenced(kcal/mol)\n"
                f"Input file: {PMF_NORMAL}\n"
                f"Rightmost reference x = {refN:.6f} Å\n"
                f"Subtracted offset = {offN:.6f} kcal/mol"
            ),
            comments="# "
        )

        np.savetxt(
            OUT_REF_CANCER,
            np.column_stack([xC, yC]),
            header=(
                "x(Å)  PMF_rightmost_referenced(kcal/mol)\n"
                f"Input file: {PMF_CANCER}\n"
                f"Rightmost reference x = {refC:.6f} Å\n"
                f"Subtracted offset = {offC:.6f} kcal/mol"
            ),
            comments="# "
        )

        print("Saved plotted normal data:", Path(OUT_REF_NORMAL).resolve())
        print("Saved plotted cancer data:", Path(OUT_REF_CANCER).resolve())

    # ---- Plot side-by-side panels ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIGSIZE, sharex=True)

    ax1.plot(xN, yN, color=COLOR_NORMAL)
    ax1.set_xlabel(XLABEL)
    ax1.set_ylabel(YLABEL)
    panel_label(ax1, "(a)")
    style_axes(ax1)
    ax1.legend(["Normal-cell-like model"], loc="best", frameon=False)

    ax2.plot(xC, yC, color=COLOR_CANCER)
    ax2.set_xlabel(XLABEL)
    panel_label(ax2, "(b)")
    style_axes(ax2)
    ax2.legend(["Cancer-associated lipid model"], loc="best", frameon=False)

    if XLIM is not None:
        ax1.set_xlim(*XLIM)
        ax2.set_xlim(*XLIM)

    if YLIM_NORMAL is not None:
        ax1.set_ylim(*YLIM_NORMAL)

    if YLIM_CANCER is not None:
        ax2.set_ylim(*YLIM_CANCER)

    fig.tight_layout()
    fig.savefig(OUT, dpi=600)
    plt.close(fig)

    print("Saved figure:", Path(OUT).resolve())
    print()
    print("Reference check:")
    print(f"Normal rightmost x = {refN:.6f} Å")
    print(f"Normal PMF at rightmost x = {yN[np.argmax(xN)]:.6f} kcal/mol")
    print(f"Cancer rightmost x = {refC:.6f} Å")
    print(f"Cancer PMF at rightmost x = {yC[np.argmax(xC)]:.6f} kcal/mol")


if __name__ == "__main__":
    main()
