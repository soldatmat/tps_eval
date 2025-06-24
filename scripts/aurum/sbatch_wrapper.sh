#!/bin/bash

# Usage: sh sbatch_wrapper.sh job_script_name [job script args...]

WRAPPER_DIR=$(dirname "$0")

SCRIPT_PATH="$WRAPPER_DIR"/jobs/"$1"
if [[ "$SCRIPT_PATH" != *.sh ]]; then
    SCRIPT_PATH="${SCRIPT_PATH}.sh"
fi

sbatch --output="$WRAPPER_DIR"/../../logs/%x.%j.out "$SCRIPT_PATH" "${@:2}"
