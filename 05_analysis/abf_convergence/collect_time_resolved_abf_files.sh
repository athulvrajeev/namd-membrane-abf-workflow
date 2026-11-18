#!/bin/bash

FIRST_WIN=2
LAST_WIN=61

BASE_OUT="PMF_convergence"
mkdir -p "$BASE_OUT"

declare -A TIME_MAP
TIME_MAP["2.0ns"]="4ns"
TIME_MAP["4.0ns"]="8ns"
TIME_MAP["6.0ns"]="12ns"
TIME_MAP["8.0ns"]="16ns"
TIME_MAP["10.0ns"]="20ns"
TIME_MAP["12.0ns"]="24ns"
TIME_MAP["14.0ns"]="28ns"
TIME_MAP["16.0ns"]="32ns"
TIME_MAP["18.0ns"]="36ns"
TIME_MAP["20.0ns"]="40ns"

for SRC_TIME in "${!TIME_MAP[@]}"; do
    DEST_TIME="${TIME_MAP[$SRC_TIME]}"
    DEST_DIR="${BASE_OUT}/PMF_${DEST_TIME}"

    mkdir -p "$DEST_DIR"

    echo "Copying $SRC_TIME -> $DEST_DIR"

    for WIN in $(seq $FIRST_WIN $LAST_WIN); do
        SRC_DIR="win${WIN}/${SRC_TIME}"

        GRAD_FILE="${SRC_DIR}/win${WIN}.grad"
        COUNT_FILE="${SRC_DIR}/win${WIN}.count"

        [ -f "$GRAD_FILE" ] && cp "$GRAD_FILE" "$DEST_DIR/" || echo "Missing: $GRAD_FILE"
        [ -f "$COUNT_FILE" ] && cp "$COUNT_FILE" "$DEST_DIR/" || echo "Missing: $COUNT_FILE"
    done
done

echo "Done. Files copied into $BASE_OUT/"
