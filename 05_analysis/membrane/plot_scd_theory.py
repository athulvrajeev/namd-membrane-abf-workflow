import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================
# Load data
# ==========================
data = pd.read_csv("SCD_graph_degrees.csv")

theta = data.iloc[:, 0]
scd = data.iloc[:, 1]

# ==========================
# House Style Settings
# ==========================
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 16,
    "axes.linewidth": 2,
    "axes.labelsize": 22,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
})

# ==========================
# Create Figure
# ==========================
fig, ax = plt.subplots(figsize=(7, 5), facecolor="white")
ax.set_facecolor("white")

# Plot line
ax.plot(theta, scd,
        color="#0b1e77",   # deep publication blue
        linewidth=3)

# Axis labels
ax.set_xlabel(r'$\theta$ (Degree)', labelpad=10)
ax.set_ylabel(r'$S_{CD}$', labelpad=10)

# Axis limits
ax.set_xlim(-20, 200)
ax.set_ylim(-0.7, 1.2)

# Ticks
ax.tick_params(direction='out', length=8, width=2)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Thicken left and bottom
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)

# Legend
ax.legend(['SCD'],
          frameon=False,
          facecolor='white',
          edgecolor='black',
          fontsize=18,
          loc='lower right')

plt.tight_layout()

# Save publication files
plt.savefig("SCD_vs_theta.pdf", dpi=600)
plt.savefig("SCD_vs_theta.png", dpi=600)

plt.show()
