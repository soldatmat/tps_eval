#!/bin/bash

############################################################
# Argument parsing                                         #
############################################################
USAGE="--cluster <cluster> --job_name <job_name> [--job_args <args...>] [--submit_args <args...>]"

Help()
{
    # Display Help
    echo "Usage: submit_job.sh $USAGE"
    echo
    echo "Arguments:"
    echo "  --cluster                   Name of the cluster"
    echo "  --job_name                  Name of the job"
    echo "  --job_args                  (optional), All following tokens (until next --option) are passed to the job script. If you need to pass arguments starting with --, enclose them in double quotes."
    echo "  --submit_args               (optional), All following tokens (until next --option) are passed to the cluster's submit command. If you need to pass arguments starting with --, enclose them in double quotes."
    echo "  -h, --help                  Show this help message and exit"
    echo
}

# Parse long options manually
JOB_ARGS=""
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --cluster)
            CLUSTER="$2"
            shift 2
            ;;
        --job_name)
            JOB_NAME="$2"
            shift 2
            ;;
        --job_args)
            shift
            # Collect all following args until the next token that starts with --
            while [[ $# -gt 0 && "$1" != --* ]]; do
                JOB_ARGS="${JOB_ARGS:+$JOB_ARGS }$1"
                shift
            done
            ;;
        --submit_args)
            shift
            # Collect all following args until the next token that starts with --
            while [[ $# -gt 0 && "$1" != --* ]]; do
                SUBMIT_ARGS="${SUBMIT_ARGS:+$SUBMIT_ARGS }$1"
                shift
            done
            ;;
        -h|--help)
            Help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            Help
            exit 1
            ;;
    esac
done

if [[ -z "$JOB_NAME" || -z "$CLUSTER" ]]; then
    echo "Both --job_name and --cluster are required."
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Remove surrounding double quotes from JOB_ARGS and SUBMIT_ARGS if present
if [[ -n "${JOB_ARGS}" ]]; then
    if [[ "${JOB_ARGS:0:1}" == "\"" && "${JOB_ARGS: -1}" == "\"" ]]; then
        JOB_ARGS="${JOB_ARGS:1:-1}"
    fi
fi
if [[ -n "${SUBMIT_ARGS}" ]]; then
    if [[ "${SUBMIT_ARGS:0:1}" == "\"" && "${SUBMIT_ARGS: -1}" == "\"" ]]; then
        SUBMIT_ARGS="${SUBMIT_ARGS:1:-1}"
    fi
fi



############################################################
# Submit job                                               #
############################################################

SCRIPTS_DIR=$(dirname "$0")

. "$SCRIPTS_DIR"/"$CLUSTER"/config.sh # Load SUBMIT_JOB variable

JOB_PATH="$SCRIPTS_DIR"/"$CLUSTER"/jobs/"$JOB_NAME"
if [[ "$JOB_PATH" != *.sh ]]; then
    JOB_PATH="${JOB_PATH}.sh"
fi
if [[ ! -f "$JOB_PATH" ]]; then
    echo "Job script $JOB_PATH does not exist."
    exit 1
fi

echo "Submitting job with command: $SUBMIT_JOB $SUBMIT_ARGS $JOB_PATH $JOB_ARGS"
$SUBMIT_JOB $SUBMIT_ARGS "$JOB_PATH" $JOB_ARGS

# The last thing printed has to be the job ID for job dependencies.
