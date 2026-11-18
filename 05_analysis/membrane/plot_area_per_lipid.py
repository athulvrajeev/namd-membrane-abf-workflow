#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ===================== USER SETTINGS =====================
NORMAL_FILE = "nonr_area_per_lipid.dat"
CANCER_FILE = "can_area_per_lipid.dat"

OUT = "Figure_APL_normal_vs_cancer.pdf"

# Locked theme colors (consistent across paper)
COLOR_NORMAL = "#2e9e59"   # green
COLOR_CANCER = "#e53935"   # red

BLOCK_NS = 10.0
ERROR = "sd"   # "sd" or "sem"

FIGSIZE = (3.6, 3.0)  # single-column width
# =========================================================


def load_two_col(path):
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def block_stats(t, y, block_ns, error="sd"):
    t0, t1 = float(t.min()), float(t.max())
    edges = np.arange(t0, t1 + block_ns, block_ns)

    t_mid, mean, err = [], [], []

    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        mask = (t >= lo) & (t < hi) if i < len(edges) - 2 else (t >= lo) & (t <= hi)

        if not np.any(mask):
            continue

        vals = y[mask]
        m = float(np.mean(vals))
        sd = float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0

        if error.lower() == "sem" and vals.size > 1:
            e = sd / np.sqrt(vals.size)
        else:
            e = sd

        t_mid.append((lo + hi) / 2.0)
        mean.append(m)
        err.append(e)

    return np.array(t_mid), np.array(mean), np.array(err)


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.tick_params(direction="out", length=5, width=1.4)


def main():
    # ---- MATCHED STYLE (same as density script) ----
    plt.rcParams.update({
        "font.family": "DejaVu Serif",
        "font.size": 10,
        "axes.labelsize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 1.5,
        "lines.linewidth": 2.2,
        "legend.fontsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    tN, yN = load_two_col(NORMAL_FILE)
    tC, yC = load_two_col(CANCER_FILE)

    tNm, mN, eN = block_stats(tN, yN, BLOCK_NS, ERROR)
    tCm, mC, eC = block_stats(tC, yC, BLOCK_NS, ERROR)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    ax.errorbar(
        tNm, mN, yerr=eN,
        fmt="-o",
        ms=3.5,
        capsize=2.5,
        color=COLOR_NORMAL,
        label="Normal membrane"
    )

    ax.errorbar(
        tCm, mC, yerr=eC,
        fmt="-o",
        ms=3.5,
        capsize=2.5,
        color=COLOR_CANCER,
        label="Cancer membrane"
    )

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(r"Area per lipid ($\AA^2$)")

    style_axes(ax)

    ax.legend(frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(OUT, dpi=600)
    plt.close(fig)

    print("Saved:", Path(OUT).resolve())


if __name__ == "__main__":
    main()
