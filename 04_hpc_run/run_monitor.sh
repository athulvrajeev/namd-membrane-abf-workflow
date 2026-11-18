#!/bin/bash
set -euo pipefail

############################################################
## USER SETTINGS
############################################################

NAMD_PATH="${NAMD_PATH:-namd3}"   # set NAMD_PATH or put namd3 on your PATH

CONFIG_BASE="$(pwd)/win"

MONITOR_SCRIPT="./monitor_copy.py"

FIRST_WINDOW=0
LAST_WINDOW=0

############################################################
## RUN WINDOWS
############################################################

for ((i=${FIRST_WINDOW}; i<=${LAST_WINDOW}; i++))
do

    WIN_DIR="${CONFIG_BASE}${i}"

    CONFIG_FILE="${WIN_DIR}/win${i}.conf"

    LOG_FILE="${WIN_DIR}/win${i}.log"

    echo "================================================="
    echo "Running ABF window ${i}"
    echo "================================================="

    ########################################################
    ## START MONITOR
    ########################################################

    python3 "${MONITOR_SCRIPT}" \
        --window "${i}" \
        --base_dir "${WIN_DIR}" \
        --output_dir_name "output" \
        --timestep_fs 1.0 &

    monitor_pid=$!

    ########################################################
    ## RUN NAMD
    ########################################################

    set +e

    "${NAMD_PATH}" \
        +ppn 18 \
        +netpoll +idlepoll \
        +devices 0 \
        "${CONFIG_FILE}" > "${LOG_FILE}" 2>&1

    namd_status=$?

    set -e

    ########################################################
    ## STOP MONITOR
    ########################################################

    kill "${monitor_pid}" 2>/dev/null || true

    wait "${monitor_pid}" 2>/dev/null || true

    ########################################################
    ## CHECK STATUS
    ########################################################

    if [ "${namd_status}" -ne 0 ]; then

        echo ""
        echo "ERROR: window ${i} failed"
        echo "Check:"
        echo "${LOG_FILE}"
        echo ""

        exit "${namd_status}"
    fi

    echo ""
    echo "Window ${i} finished successfully"
    echo ""

done

echo "=========================================="
echo "ALL WINDOWS COMPLETED"
echo "=========================================="
