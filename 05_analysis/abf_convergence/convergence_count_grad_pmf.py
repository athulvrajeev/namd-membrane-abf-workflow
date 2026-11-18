#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================

SYSTEMS = {
    "bare_normal": {
        "abf_dir": "data/ABF/bare_normal",
        "pmf_dir": "data/ABF/bare_normal/PMF_convergence/pmf_convergence",
        "pmf_prefix": "merge_nor_gold",
    },
    "bare_cancer": {
        "abf_dir": "data/ABF/bare_cancer",
        "pmf_dir": "data/ABF/bare_cancer/PMF_convergence/pmf_convergence",
        "pmf_prefix": "merge_can_gold",
    },
    "meth_normal": {
        "abf_dir": "data/ABF/Meth_normal",
        "pmf_dir": "data/ABF/Meth_normal/PMF_convergence/pmf_convergence",
        "pmf_prefix": "merge_meth_normal",
    },
    "meth_cancer": {
        "abf_dir": "data/ABF/Meth_cancer",
        "pmf_dir": "data/ABF/Meth_cancer/PMF_convergence/pmf_convergence",
        "pmf_prefix": "merge_meth_cancer",
    },
    "oct_normal": {
        "abf_dir": "data/ABF/Oct_normal",
        "pmf_dir": "data/ABF/Oct_normal/PMF_convergence/pmf_convergence",
        "pmf_prefix": "merge_oct_normal",
    },
    "oct_cancer": {
        "abf_dir": "data/ABF/oct_cancer",
        "pmf_dir": "data/ABF/oct_cancer/PMF_convergence/pmf_convergence",
        "pmf_prefix": "merge_oct_can",
    },
}

FIRST_WIN = 0
LAST_WIN = 69

TIME_MAP = [
    ("2.0ns", 4),
    ("4.0ns", 8),
    ("6.0ns", 12),
    ("8.0ns", 16),
    ("10.0ns", 20),
    ("12.0ns", 24),
    ("14.0ns", 28),
    ("16.0ns", 32),
    ("18.0ns", 36),
    ("20.0ns", 40),
]

# For count plots, fewer late-time curves are usually clearer
COUNT_TIME_MAP = [
    ("14.0ns", 28),
    ("16.0ns", 32),
    ("18.0ns", 36),
    ("20.0ns", 40),
]

OUTPUT_DIR = "SI_three_panel_convergence"

REFERENCE_MODE = "first"
# Options:
# "first" -> set first PMF point to zero
# "last"  -> set last PMF point to zero
# float   -> e.g., 30.0 sets PMF at coordinate closest to 30 Å to zero

FIGSIZE = (7.2, 9.0)

LINEWIDTH = 1.7
MARKERSIZE = 2.8
MARK_EVERY_GRAD = 8
MARK_EVERY_COUNT = 8
MARK_EVERY_PMF = 25

XLABEL = r"Reaction coordinate ($\mathrm{\AA}$)"
YLABEL_COUNT = r"Cumulative count ($\times 10^6$)"
YLABEL_GRAD = r"Mean force (kcal mol$^{-1}$ $\mathrm{\AA}^{-1}$)"
YLABEL_PMF = r"PMF (kcal mol$^{-1}$)"

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
        "axes.labelsize": 11,
        "axes.titlesize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 1.4,

        "legend.fontsize": 10,
        "lines.linewidth": LINEWIDTH,
    })


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.minorticks_off()
    ax.tick_params(direction="out", length=5, width=1.2)
    ax.axvline(
        0,
        color="0.55",
        linestyle="--",
        linewidth=1.0,
        alpha=0.85,
    )


def load_merged_window_data(
    base_folder: str,
    folder_name: str,
    ext: str,
):
    all_data = []

    for win in range(FIRST_WIN, LAST_WIN + 1):
        fpath = (
            Path(base_folder)
            / f"win{win}"
            / folder_name
            / f"win{win}.{ext}"
        )

        if not fpath.exists():
            continue

        try:
            data = np.loadtxt(fpath)

            if data.ndim == 2 and data.shape[1] >= 2:
                all_data.append(data[:, :2])

        except Exception as e:
            print(f"    [Warning] Could not read {fpath}: {e}")

    if not all_data:
        return None

    merged = np.vstack(all_data)
    merged = merged[np.argsort(merged[:, 0])]

    return merged


def read_pmf(fname: Path):
    data = []

    with open(fname, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 2:
                continue

            try:
                data.append([
                    float(parts[0]),
                    float(parts[1]),
                ])
            except ValueError:
                continue

    if not data:
        return None

    arr = np.asarray(data, dtype=float)
    arr = arr[np.argsort(arr[:, 0])]

    return arr[:, 0], arr[:, 1]


def apply_reference(x, y, mode):
    y = y.copy()

    if mode == "first":
        return y - y[0]

    if mode == "last":
        return y - y[-1]

    ref_x = float(mode)
    idx = np.argmin(np.abs(x - ref_x))

    return y - y[idx]


def calculate_late_pmf_drift(curves):
    """
    Compare 32 ns and 40 ns if both exist.
    Otherwise compare the final two available curves.
    """
    if len(curves) < 2:
        return None

    by_time = {
        t: (x, y)
        for t, x, y in curves
    }

    if 32 in by_time and 40 in by_time:
        t1, t2 = 32, 40
        x1, y1 = by_time[t1]
        x2, y2 = by_time[t2]

    else:
        t1, x1, y1 = curves[-2]
        t2, x2, y2 = curves[-1]

    if not np.allclose(x1, x2):
        print(
            "    [Warning] PMF grids do not match; "
            "interpolating for drift calculation."
        )

        y2_interp = np.interp(x1, x2, y2)
        diff = y2_interp - y1

    else:
        diff = y2 - y1

    max_abs = float(np.max(np.abs(diff)))
    rms = float(np.sqrt(np.mean(diff**2)))

    return t1, t2, max_abs, rms


def plot_system(label: str, info: dict):
    abf_dir = info["abf_dir"]
    pmf_dir = Path(info["pmf_dir"])
    pmf_prefix = info["pmf_prefix"]

    print(f"\n[{label}]")
    print(f"  ABF dir: {abf_dir}")
    print(f"  PMF dir: {pmf_dir}")

    fig, axes = plt.subplots(
        3,
        1,
        figsize=FIGSIZE,
        sharex=False,
        gridspec_kw={
            "height_ratios": [1, 1.15, 1.15]
        },
    )

    ax_count, ax_grad, ax_pmf = axes

    colors = [
        "#4d4d4d",
        "#1565c0",
        "#2e9e59",
        "#a970d6",
        "#e67e22",
        "#8e44ad",
        "#00a6a6",
        "#d4a017",
        "#7f8c8d",
        "#e53935",
    ]

    markers = [
        "o",
        "^",
        "s",
        "D",
        "p",
        "h",
        "v",
        "<",
        ">",
        "*",
    ]

    # --------------------------------------------------------
    # Panel A: Count convergence
    # --------------------------------------------------------
    found_count = False

    start_idx = len(colors) - len(COUNT_TIME_MAP)

    for i, (folder_name, real_ns) in enumerate(COUNT_TIME_MAP):
        merged = load_merged_window_data(
            abf_dir,
            folder_name,
            "count",
        )

        if merged is None:
            print(
                f"    [Info] Missing count data for {real_ns} ns"
            )
            continue

        found_count = True
        ci = start_idx + i

        ax_count.plot(
            merged[:, 0],
            merged[:, 1] / 1e6,
            color=colors[ci],
            marker=markers[ci],
            markersize=MARKERSIZE,
            markevery=MARK_EVERY_COUNT,
            linewidth=LINEWIDTH,
            label=(
                f"{real_ns} ns"
                + (
                    " final"
                    if real_ns == COUNT_TIME_MAP[-1][1]
                    else ""
                )
            ),
        )

    ax_count.axhline(
        1000 / 1e6,
        color="0.45",
        linestyle=":",
        linewidth=1.0,
        alpha=0.9,
        label="fullSamples",
    )

    ax_count.set_ylabel(YLABEL_COUNT)

    ax_count.text(
        0.02,
        0.90,
        "(a)",
        transform=ax_count.transAxes,
        fontweight="bold",
    )

    style_axes(ax_count)

    if found_count:
        ax_count.legend(
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.28),
            ncol=5,
            columnspacing=0.8,
            handlelength=1.8,
        )

    # --------------------------------------------------------
    # Panel B: Grad convergence
    # --------------------------------------------------------
    found_grad = False

    for i, (folder_name, real_ns) in enumerate(TIME_MAP):
        merged = load_merged_window_data(
            abf_dir,
            folder_name,
            "grad",
        )

        if merged is None:
            print(
                f"    [Info] Missing grad data for {real_ns} ns"
            )
            continue

        found_grad = True

        ax_grad.plot(
            merged[:, 0],
            merged[:, 1],
            color=colors[i],
            marker=markers[i],
            markersize=MARKERSIZE,
            markevery=MARK_EVERY_GRAD,
            linewidth=LINEWIDTH,
            label=(
                f"{real_ns} ns"
                + (
                    " final"
                    if real_ns == TIME_MAP[-1][1]
                    else ""
                )
            ),
        )

    ax_grad.axhline(
        0,
        color="0.45",
        linestyle=":",
        linewidth=1.0,
        alpha=0.9,
    )

    ax_grad.set_ylabel(YLABEL_GRAD)

    ax_grad.text(
        0.02,
        0.90,
        "(b)",
        transform=ax_grad.transAxes,
        fontweight="bold",
    )

    style_axes(ax_grad)

    if found_grad:
        ax_grad.legend(
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.34),
            ncol=5,
            columnspacing=0.75,
            handlelength=1.7,
        )

    # --------------------------------------------------------
    # Panel C: PMF convergence
    # --------------------------------------------------------
    pmf_curves = []

    for i, (_, real_ns) in enumerate(TIME_MAP):
        fname = pmf_dir / f"{pmf_prefix}_{real_ns}ns.pmf"

        if not fname.exists():
            print(f"    [Info] Missing PMF file: {fname}")
            continue

        pmf_data = read_pmf(fname)

        if pmf_data is None:
            print(f"    [Info] Could not read PMF file: {fname}")
            continue

        x, y = pmf_data
        y_ref = apply_reference(
            x,
            y,
            REFERENCE_MODE,
        )

        ax_pmf.plot(
            x,
            y_ref,
            color=colors[i],
            marker=markers[i],
            markersize=MARKERSIZE,
            markevery=MARK_EVERY_PMF,
            linewidth=LINEWIDTH,
            label=(
                f"{real_ns} ns"
                + (
                    " final"
                    if real_ns == TIME_MAP[-1][1]
                    else ""
                )
            ),
        )

        pmf_curves.append((
            real_ns,
            x,
            y_ref,
        ))

    ax_pmf.axhline(
        0,
        color="0.45",
        linestyle=":",
        linewidth=1.0,
        alpha=0.9,
    )

    ax_pmf.set_xlabel(XLABEL)
    ax_pmf.set_ylabel(YLABEL_PMF)

    ax_pmf.text(
        0.02,
        0.90,
        "(c)",
        transform=ax_pmf.transAxes,
        fontweight="bold",
    )

    style_axes(ax_pmf)

    if pmf_curves:
        ax_pmf.legend(
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.34),
            ncol=5,
            columnspacing=0.75,
            handlelength=1.7,
        )

    # --------------------------------------------------------
    # Final figure layout — no overall title
    # --------------------------------------------------------
    fig.tight_layout(h_pad=2.2)

    out_base = (
        Path(OUTPUT_DIR)
        / f"SI_three_panel_convergence_{label}"
    )

    fig.savefig(
        f"{out_base}.pdf",
        bbox_inches="tight",
    )

    fig.savefig(
        f"{out_base}.png",
        dpi=600,
        bbox_inches="tight",
    )

    plt.close(fig)

    drift = calculate_late_pmf_drift(pmf_curves)

    if drift is not None:
        t1, t2, max_abs, rms = drift

        print(
            f"  PMF drift {t1} vs {t2} ns: "
            f"Max |ΔPMF| = {max_abs:.3f}, "
            f"RMS = {rms:.3f}"
        )

        with open(
            f"{out_base}_metrics.txt",
            "w",
        ) as f:
            f.write(f"System: {label}\n")
            f.write(
                f"Compared PMF profiles: "
                f"{t1} ns vs {t2} ns\n"
            )
            f.write(
                f"Max |Delta PMF|: "
                f"{max_abs:.6f} kcal/mol\n"
            )
            f.write(
                f"RMS Delta PMF: "
                f"{rms:.6f} kcal/mol\n"
            )

    print(f"  Saved: {out_base}.pdf")
    print(f"  Saved: {out_base}.png")


def main():
    set_pccp_style()
    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    print("=" * 70)
    print("Generating three-panel ABF convergence figures")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 70)

    for label, info in SYSTEMS.items():
        plot_system(
            label,
            info,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()

