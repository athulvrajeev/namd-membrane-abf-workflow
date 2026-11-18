#!/usr/bin/env python3
## ============================================================
## run_all_convergence.py
## Batch runner — generates SI convergence figures for all
## 6 systems. Give the full path to each system's base folder
## in the USER SETTINGS section below.
##
## Expected structure inside each base folder (uniform):
##   <base_folder>/win0/2.0ns/win0.grad
##   <base_folder>/win0/2.0ns/win0.count
##   <base_folder>/win1/2.0ns/win1.grad
##   ...
##   <base_folder>/win69/20.0ns/win69.grad
##
## Usage: python run_all_convergence.py
## ============================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
import os

## ============================================================
## USER SETTINGS — set full paths for each system here
## ============================================================

## Format: "label": ("/full/path/to/base/folder", "Display Name")
## Replace the paths below with your actual folder paths.

SYSTEMS = {
    "bare_normal": (
        "/path/to/bare_gold/normal_membrane",
        "Bare AuNP — Normal-cell-like model"
    ),
    "bare_cancer": (
        "/path/to/bare_gold/cancer_membrane",
        "Bare AuNP — Cancer-associated lipid model"
    ),
    "meth_normal": (
        "/path/to/methane_thiol/normal_membrane",
        "Methane Thiol AuNP — Normal-cell-like model"
    ),
    "meth_cancer": (
        "/path/to/methane_thiol/cancer_membrane",
        "Methane Thiol AuNP — Cancer-associated lipid model"
    ),
    "oct_normal": (
        "/path/to/octane_thiol/normal_membrane",
        "Octane Thiol AuNP — Normal-cell-like model"
    ),
    "oct_cancer": (
        "/path/to/octane_thiol/cancer_membrane",
        "Octane Thiol AuNP — Cancer-associated lipid model"
    ),
}

FIRST_WIN = 0
LAST_WIN  = 69

FOLDER_TIME_PAIRS = [
    ("2.0ns",  4),  ("4.0ns",  8),  ("6.0ns",  12),
    ("8.0ns",  16), ("10.0ns", 20), ("12.0ns", 24),
    ("14.0ns", 28), ("16.0ns", 32), ("18.0ns", 36),
    ("20.0ns", 40),
]
COUNT_TIME_PAIRS = [("16.0ns", 32), ("18.0ns", 36), ("20.0ns", 40)]

OUTPUT_DIR = "SI_convergence_figures"

## ============================================================

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11,
    "axes.linewidth": 1.2,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "in", "ytick.direction": "in",
})

os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_merged(base_folder, win_range, folder_name, ext):
    all_data = []
    for win in range(win_range[0], win_range[1] + 1):
        fpath = Path(base_folder) / f"win{win}" / folder_name / f"win{win}.{ext}"
        if not fpath.exists():
            continue
        try:
            data = np.loadtxt(fpath)
            if data.ndim == 2 and data.shape[1] >= 2:
                all_data.append(data[:, :2])
        except Exception:
            continue
    if not all_data:
        return None
    merged = np.vstack(all_data)
    return merged[np.argsort(merged[:, 0])]

def make_grad_figure(label, base_folder, display_name):
    fig, ax = plt.subplots(figsize=(10, 6))
    cmap   = cm.plasma
    ntimes = len(FOLDER_TIME_PAIRS)
    colors = [cmap(0.15 + 0.75 * i / max(1, ntimes-1)) for i in range(ntimes)]

    found_any = False
    for i, (folder_name, real_ns) in enumerate(FOLDER_TIME_PAIRS):
        merged = load_merged(base_folder, (FIRST_WIN, LAST_WIN), folder_name, "grad")
        if merged is None:
            continue
        found_any = True
        lw    = 2.2 if i == ntimes-1 else 1.0
        alpha = 1.0 if i == ntimes-1 else 0.75
        zord  = 5   if i == ntimes-1 else 2
        lbl   = f"{real_ns} ns" + (" (final)" if i == ntimes-1 else "")
        ax.plot(merged[:, 0], merged[:, 1],
                color=colors[i], linewidth=lw,
                alpha=alpha, zorder=zord, label=lbl)

    if not found_any:
        print(f"  [Warning] No .grad files found for {label} in {base_folder}")
        plt.close()
        return

    ax.axhline(0, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', linewidth=1.0,
               alpha=0.6, label='Membrane centre')
    ax.set_xlabel("Reaction coordinate (Å)", fontsize=12)
    ax.set_ylabel("Mean force (kcal mol$^{-1}$ Å$^{-1}$)", fontsize=12)
    ax.set_title(f"ABF Gradient Convergence — {display_name}",
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, ncol=2, framealpha=0.8, edgecolor='lightgray')
    ax.grid(True, alpha=0.2, linestyle=':')
    plt.tight_layout()

    out = Path(OUTPUT_DIR) / f"SI_grad_convergence_{label}"
    fig.savefig(f"{out}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: SI_grad_convergence_{label}.pdf")

def make_count_figure(label, base_folder, display_name):
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(f"ABF Sampling — {display_name}",
                 fontsize=11, fontweight='bold')

    cmap2   = cm.viridis
    ntimes2 = len(COUNT_TIME_PAIRS)
    colors2 = [cmap2(0.2 + 0.6*i / max(1, ntimes2-1)) for i in range(ntimes2)]

    found_any = False
    prev_merged = None
    prev_ns     = None
    for i, (folder_name, real_ns) in enumerate(COUNT_TIME_PAIRS):
        merged = load_merged(base_folder, (FIRST_WIN, LAST_WIN), folder_name, "count")
        if merged is None:
            print(f"    [Info] No count data at {real_ns} ns for {label}")
            continue
        found_any = True
        lw  = 2.0 if i == ntimes2-1 else 1.2
        ax_a.plot(merged[:, 0], merged[:, 1] / 1e6,
                  color=colors2[i], linewidth=lw, label=f"{real_ns} ns")

        ## Incremental panel — interpolate onto common grid to handle
        ## arrays of different lengths across time points
        if prev_merged is not None:
            xi_curr = merged[:, 0]
            xi_prev = prev_merged[:, 0]
            xi_min  = max(xi_curr.min(), xi_prev.min())
            xi_max  = min(xi_curr.max(), xi_prev.max())
            if xi_max > xi_min:
                xi_common = np.linspace(xi_min, xi_max, 500)
                cnt_curr  = np.interp(xi_common, xi_curr, merged[:, 1])
                cnt_prev  = np.interp(xi_common, xi_prev, prev_merged[:, 1])
                diff      = cnt_curr - cnt_prev
                ax_b.plot(xi_common, diff / 1e3,
                          color=colors2[i], linewidth=lw,
                          label=f"{prev_ns}→{real_ns} ns")

        prev_merged = merged
        prev_ns     = real_ns

    if not found_any:
        print(f"  [Warning] No .count files found for {label} in {base_folder}")
        plt.close()
        return

    ax_a.axhline(1000/1e6, color='red', linestyle='--', linewidth=1.2,
                 label='fullSamples threshold (1000)', zorder=5)
    for ax in [ax_a, ax_b]:
        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.grid(True, alpha=0.2, linestyle=':')
    ax_b.axhline(0, color='gray', linestyle=':', linewidth=0.8)

    ax_a.set_ylabel("Cumulative count (×10⁶)", fontsize=11)
    ax_a.set_title("Cumulative sampling per bin", fontsize=10)
    ax_a.legend(fontsize=9, framealpha=0.8)

    ax_b.set_xlabel("Reaction coordinate (Å)", fontsize=12)
    ax_b.set_ylabel("New samples per interval (×10³)", fontsize=11)
    ax_b.set_title("Incremental sampling between time points", fontsize=10)
    ax_b.legend(fontsize=9, framealpha=0.8)

    plt.tight_layout()
    out = Path(OUTPUT_DIR) / f"SI_count_convergence_{label}"
    fig.savefig(f"{out}.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{out}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: SI_count_convergence_{label}.pdf")

## ============================================================
## MAIN
## ============================================================

print("=" * 60)
print("Generating SI convergence figures for all systems")
print(f"Output directory: {OUTPUT_DIR}/")
print("=" * 60)

for label, (base_folder, display_name) in SYSTEMS.items():
    print(f"\n[{label}]")
    if not Path(base_folder).is_dir():
        print(f"  [Skipping] Folder not found: {base_folder}")
        continue
    make_grad_figure(label,  base_folder, display_name)
    make_count_figure(label, base_folder, display_name)

print("\n" + "=" * 60)
print(f"Done. All figures saved to: {OUTPUT_DIR}/")
print("Files generated per system:")
print("  SI_grad_convergence_<label>.pdf")
print("  SI_count_convergence_<label>.pdf")
print("=" * 60)
