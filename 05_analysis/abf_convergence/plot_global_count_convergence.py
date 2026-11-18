#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

FIRST_WIN = 0
LAST_WIN = 69

# Plot only later-time count
folder_real_time_pairs = [
#    ("10.0ns", 20),
  #  ("12.0ns", 24),
  #  ("14.0ns", 28),
    ("16.0ns", 32),
    ("18.0ns", 36),
    ("20.0ns", 40),
]

plt.figure(figsize=(12, 8))

for folder_name, real_time_ns in folder_real_time_pairs:

    all_data = []

    for win in range(FIRST_WIN, LAST_WIN + 1):

        count_file = Path(f"win{win}") / folder_name / f"win{win}.count"

        if not count_file.exists():
            print(f"Missing: {count_file}")
            continue

        data = np.loadtxt(count_file)

        if data.ndim == 2 and data.shape[1] >= 2:
            all_data.append(data[:, :2])

    if not all_data:
        print(f"No data for {real_time_ns} ns")
        continue

    merged = np.vstack(all_data)
    merged = merged[np.argsort(merged[:, 0])]

    plt.plot(merged[:, 0], merged[:, 1], label=f"{real_time_ns} ns", linewidth=1.4)

plt.xlabel("Reaction coordinate (Å)")
plt.ylabel("Count")
plt.title("Global ABF Count Convergence: 20–40 ns")
plt.legend(fontsize=9)
plt.tight_layout()

plt.savefig("global_count_convergence_20to40ns.png", dpi=300)
print("Saved global_count_convergence_20to40ns.png")
