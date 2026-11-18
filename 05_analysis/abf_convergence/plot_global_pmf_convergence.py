#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

FIRST_WIN = 0
LAST_WIN = 69

# Because timestep is 2 fs:
# folder 2.0ns = real 4 ns, folder 4.0ns = real 8 ns, etc.
folder_real_time_pairs = [
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

plt.figure(figsize=(12, 8))

for folder_name, real_time_ns in folder_real_time_pairs:

    all_data = []

    for win in range(FIRST_WIN, LAST_WIN + 1):

        pmf_file = Path(f"win{win}") / folder_name / f"win{win}.pmf"

        if not pmf_file.exists():
            print(f"Missing: {pmf_file}")
            continue

        data = np.loadtxt(pmf_file)

        if data.ndim == 2 and data.shape[1] >= 2:
            all_data.append(data[:, :2])

    if not all_data:
        print(f"No data for {real_time_ns} ns")
        continue

    merged = np.vstack(all_data)

    # Sort by reaction coordinate
    merged = merged[np.argsort(merged[:, 0])]

    plt.plot(
        merged[:, 0],
        merged[:, 1],
        label=f"{real_time_ns} ns",
        linewidth=1.2
    )

plt.xlabel("Reaction coordinate (Å)")
plt.ylabel("Potential of Mean force (kcal/mol)")
plt.title("Global ABF pmf Convergence")
plt.legend(fontsize=8)
plt.tight_layout()

plt.savefig("global_pmf_convergence_4to40ns.png", dpi=300)
print("Saved global_pmf_convergence_4to40ns.png")
