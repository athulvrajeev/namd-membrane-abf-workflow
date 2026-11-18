#!/usr/bin/env python3

import time
import shutil
import argparse
from pathlib import Path

############################################################
## ARGUMENTS
############################################################

parser = argparse.ArgumentParser()

parser.add_argument("--window", type=int, required=True)

parser.add_argument("--base_dir", required=True)

parser.add_argument("--output_dir_name", default="output")

parser.add_argument("--timestep_fs", type=float, default=1.0)

args = parser.parse_args()

############################################################
## SETTINGS
############################################################

win = args.window

base_dir = Path(args.base_dir)

output_dir = base_dir / args.output_dir_name

############################################################
## SAVE EVERY 1 ns
############################################################

interval_steps = 1000000

############################################################
## FILES TO COPY
############################################################

source_files = [
    f"win{win}.colvars.state",
    f"win{win}.colvars.traj",
    f"win{win}.count",
    f"win{win}.grad",
    f"win{win}.xst",
    f"win{win}.pmf"
]

xst_file = output_dir / f"win{win}.xst"

next_target = interval_steps

############################################################
## GET LAST STEP
############################################################

def get_last_step():

    try:

        with open(xst_file) as f:

            lines = f.readlines()

            for line in reversed(lines):

                if not line.startswith("#"):

                    return int(line.split()[0])

    except:
        return None

############################################################
## MAIN LOOP
############################################################

while True:

    ########################################################
    ## WAIT FOR XST
    ########################################################

    if not xst_file.exists():

        time.sleep(10)

        continue

    ########################################################
    ## GET CURRENT STEP
    ########################################################

    step = get_last_step()

    if step is None:

        time.sleep(10)

        continue

    ########################################################
    ## COPY FILES EVERY 1 ns
    ########################################################

    if step >= next_target:

        ####################################################
        ## STEP -> TIME
        ####################################################

        time_ps = step * args.timestep_fs / 1000.0

        time_ns = time_ps / 1000.0

        folder_name = f"{time_ns:.1f}ns"

        target_folder = base_dir / folder_name

        target_folder.mkdir(exist_ok=True)

        print(f"[win{win}] Saving files at {time_ns:.1f} ns")

        ####################################################
        ## COPY FILES
        ####################################################

        for file in source_files:

            src = output_dir / file

            if src.exists():

                shutil.copy2(src, target_folder / file)

        next_target += interval_steps

    ########################################################
    ## WAIT
    ########################################################

    time.sleep(10)
