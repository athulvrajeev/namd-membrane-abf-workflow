#!/usr/bin/env python3

from pathlib import Path
import numpy as np

# ================= USER SETTINGS =================
PMF_FILES = [
    "merge_oct_normal.pmf",
    "merge_oct_can.pmf",
]

OUTPUT_SUFFIX = "_right0_normalized.pmf"
# =================================================


def load_numeric_pmf(filename):
    """
    Reads PMF file and extracts numeric rows.
    First column = reaction coordinate
    Second column = PMF
    Extra numeric columns, if present, are preserved.
    Comment/header lines are ignored.
    """
    rows = []

    with open(filename, "r") as f:
        for line in f:
            s = line.strip()

            if not s:
                continue

            if s.startswith(("#", "@", ";")):
                continue

            parts = s.replace(",", " ").split()

            values = []
            for p in parts:
                try:
                    values.append(float(p))
                except ValueError:
                    pass

            if len(values) >= 2:
                rows.append(values)

    if not rows:
        raise ValueError(f"No numeric PMF data found in {filename}")

    return np.array(rows, dtype=float)


def normalize_rightmost(filename):
    filename = Path(filename)

    data = load_numeric_pmf(filename)

    x = data[:, 0]
    pmf = data[:, 1]

    # Find rightmost x-coordinate
    right_index = np.argmax(x)
    right_x = x[right_index]
    right_pmf = pmf[right_index]

    # Normalize PMF so that rightmost point becomes zero
    data[:, 1] = pmf - right_pmf

    output_file = filename.with_name(filename.stem + OUTPUT_SUFFIX)

    header = (
        f"Right-reference normalized PMF\n"
        f"Input file: {filename.name}\n"
        f"Reference point: x = {right_x:.6f} Å\n"
        f"Original PMF at reference point = {right_pmf:.6f} kcal/mol\n"
        f"Normalized PMF = PMF - PMF(rightmost x)\n"
        f"Columns: x  PMF_normalized"
    )

    np.savetxt(
        output_file,
        data,
        fmt="%.8f",
        header=header,
        comments="# "
    )

    print(f"Processed: {filename}")
    print(f"  Rightmost x = {right_x:.6f} Å")
    print(f"  Original PMF at rightmost x = {right_pmf:.6f} kcal/mol")
    print(f"  Output: {output_file}")
    print(f"  Check normalized PMF at rightmost x = {data[right_index, 1]:.6f} kcal/mol")
    print()


def main():
    for pmf_file in PMF_FILES:
        normalize_rightmost(pmf_file)


if __name__ == "__main__":
    main()
