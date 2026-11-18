#!/bin/bash

# Set input files
LOG1="../smd.log"
LOG2="../smd_1.log"
MERGED_LOG="merged_smd.log"

# Merge the two logs into one file
cat "$LOG1" "$LOG2" > "$MERGED_LOG"

echo "Merged $LOG1 and $LOG2 into $MERGED_LOG"

