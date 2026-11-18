#!/bin/bash
set -euo pipefail

BASE_DIR="$(pwd)"   # run from the system directory
INPUT_DIR="${BASE_DIR}/input"

for ((i=2; i<=61; i++))
do
    WIN_DIR="${BASE_DIR}/win${i}"

    cp "${INPUT_DIR}/Can_gold-${i}.restart.coor" "${WIN_DIR}/"
    cp "${INPUT_DIR}/Can_gold-${i}.restart.xsc"  "${WIN_DIR}/"

    echo "Copied restart files to win${i}"
done

echo "Done."
