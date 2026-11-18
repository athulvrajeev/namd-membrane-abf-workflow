#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.ndimage import gaussian_filter1d

# --------------------------------------------------
# Edit file names here
# --------------------------------------------------
systems = {
    "Bare AuNP": {
        "normal_water": "water_core_count_bare_normal.dat",
        "cancer_water": "water_core_and_local_count_bare_cancer.dat",
        "normal_thick": "local_thickness_bare_normal.dat",
        "cancer_thick": "local_thickness_bare_cancer.dat",
    },
    "Methane thiol AuNP": {
        "normal_water": "water_core_count_meth_normal.dat",
        "cancer_water": "water_core_count_meth_cancer.dat",
        "normal_thick": "local_thickness_meth_normal.dat",
        "cancer_thick": "local_thickness_meth_cancer.dat",
    },
    "Octane thiol AuNP": {
        "normal_water": "water_core_count_oct_normal.dat",
        "cancer_water": "water_core_count_oct_cancer.dat",
        "normal_thick": "local_thickness_oct_normal.dat",
        "cancer_thick": "local_thickness_oct_cancer.dat",
    },
}

# --------------------------------------------------
# Plot colours
# --------------------------------------------------
NORMAL_COLOUR = "#1565c0"
CANCER_COLOUR = "#e53935"


def load_dat(filename):
    filename = Path(filename)

    if not filename.exists():
        print(f"Missing file: {filename}")
        return None, None

    data = np.genfromtxt(filename, comments="#", dtype=float)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    z = data[:, 1]
    y = data[:, 2]

    mask = np.isfinite(z) & np.isfinite(y)

    return z[mask], y[mask]


def bin_average(z, y, bin_width=1.0):
    z_min = np.floor(np.nanmin(z))
    z_max = np.ceil(np.nanmax(z))

    bins = np.arange(z_min, z_max + bin_width, bin_width)
    centers = 0.5 * (bins[:-1] + bins[1:])

    mean_y = []
    sem_y = []

    for i in range(len(bins) - 1):
        mask = (z >= bins[i]) & (z < bins[i + 1])
        values = y[mask]

        if len(values) > 1:
            mean_y.append(np.mean(values))
            sem_y.append(np.std(values, ddof=1) / np.sqrt(len(values)))
        elif len(values) == 1:
            mean_y.append(values[0])
            sem_y.append(0.0)
        else:
            mean_y.append(np.nan)
            sem_y.append(np.nan)

    return centers, np.array(mean_y), np.array(sem_y)


def smooth_nan_array(y, sigma=1.0):
    """
    Smooth array containing NaN values without spreading NaNs.
    """
    y = np.asarray(y, dtype=float)

    valid = np.isfinite(y)
    if np.sum(valid) < 3:
        return y

    y_filled = np.where(valid, y, 0.0)
    weights = valid.astype(float)

    smooth_y = gaussian_filter1d(y_filled, sigma=sigma)
    smooth_w = gaussian_filter1d(weights, sigma=sigma)

    with np.errstate(invalid="ignore", divide="ignore"):
        result = smooth_y / smooth_w

    result[smooth_w == 0] = np.nan
    return result


def style_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=4, width=1)
    ax.grid(False)


def plot_metric(
    metric_type,
    ylabel,
    output_name,
    bin_width=1.0,
    smooth=True,
    smooth_sigma=1.0,
):
    fig, axes = plt.subplots(3, 1, figsize=(7, 12), sharex=False, sharey=False)

    panel_labels = ["(a)", "(b)", "(c)"]

    for ax, (title, files), panel in zip(axes, systems.items(), panel_labels):

        if metric_type == "water":
            normal_file = files["normal_water"]
            cancer_file = files["cancer_water"]
        elif metric_type == "thickness":
            normal_file = files["normal_thick"]
            cancer_file = files["cancer_thick"]
        else:
            raise ValueError("metric_type must be 'water' or 'thickness'")

        zn, yn = load_dat(normal_file)
        zc, yc = load_dat(cancer_file)

        if zn is not None:
            zn_bin, yn_mean, yn_sem = bin_average(zn, yn, bin_width=bin_width)

            if smooth:
                yn_plot = smooth_nan_array(yn_mean, sigma=smooth_sigma)
            else:
                yn_plot = yn_mean

            ax.plot(
                zn_bin,
                yn_plot,
                label="Normal membrane model",
                linewidth=2.2,
                color=NORMAL_COLOUR,
            )

            ax.fill_between(
                zn_bin,
                yn_plot - yn_sem,
                yn_plot + yn_sem,
                color=NORMAL_COLOUR,
                alpha=0.18,
                linewidth=0,
            )

        if zc is not None:
            zc_bin, yc_mean, yc_sem = bin_average(zc, yc, bin_width=bin_width)

            if smooth:
                yc_plot = smooth_nan_array(yc_mean, sigma=smooth_sigma)
            else:
                yc_plot = yc_mean

            ax.plot(
                zc_bin,
                yc_plot,
                label="Cancer-associated lipid model",
                linewidth=2.2,
                color=CANCER_COLOUR,
            )

            ax.fill_between(
                zc_bin,
                yc_plot - yc_sem,
                yc_plot + yc_sem,
                color=CANCER_COLOUR,
                alpha=0.18,
                linewidth=0,
            )

        ax.text(
            0.02,
            0.88,
            panel,
            transform=ax.transAxes,
            fontsize=14,
            fontweight="bold",
        )

        ax.set_title(title, fontsize=14, pad=8)
        ax.set_xlabel("NP z-position relative to membrane COM (Å)", fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.legend(frameon=False, fontsize=11)
        style_axis(ax)

    plt.tight_layout()
    plt.savefig(output_name, dpi=600, bbox_inches="tight")
    plt.savefig(output_name.replace(".png", ".pdf"), bbox_inches="tight")
    plt.show()


# --------------------------------------------------
# Generate water-count figure
# --------------------------------------------------
plot_metric(
    metric_type="water",
    ylabel="Water count in membrane-core cylinder",
    output_name="Fig_water_core_count_binned.png",
    bin_width=1.0,
    smooth=True,
    smooth_sigma=1.0,
)


# --------------------------------------------------
# Generate local-thickness figure
# --------------------------------------------------
plot_metric(
    metric_type="thickness",
    ylabel="Local membrane thickness (Å)",
    output_name="Fig_local_membrane_thickness_binned.png",
    bin_width=1.0,
    smooth=True,
    smooth_sigma=1.0,
)
