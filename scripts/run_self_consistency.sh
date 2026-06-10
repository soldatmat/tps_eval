#!/bin/bash

USAGE="--structs_dir <structs_dir> [--save_path <save_path>] [--num_seqs <n>] [--sampling_temp <t>] [--model_name <name>] [--seed <n>] [--ids <id...>] [--limit <n>] [--device <cuda|cpu>]"

Help()
{
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --structs_dir     Directory of structures (AF3 af_output or flat .pdb/.cif); stem = ID (required)"
    echo "  --save_path       Output CSV path (optional; default <structs_dir>_self_consistency.csv)"
    echo "  --num_seqs        ProteinMPNN sequences sampled per backbone (optional; default 8)"
    echo "  --sampling_temp   ProteinMPNN sampling temperature (optional; default 0.1)"
    echo "  --model_name      ProteinMPNN model name (optional; default v_48_020)"
    echo "  --seed            Random seed (optional; default 0 = random)"
    echo "  --ids             Restrict to these structure IDs (optional; validate on 1-2 first)"
    echo "  --limit           Score only the first N structures (optional; cheap validation)"
    echo "  --device          Torch device cuda/cpu for ESMFold (optional)"
    echo "  -h, --help        Show this help message and exit"
    echo
}

ids=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --structs_dir)
            structs_dir="$2"
            shift 2
            ;;
        --save_path)
            save_path="$2"
            shift 2
            ;;
        --num_seqs)
            num_seqs="$2"
            shift 2
            ;;
        --sampling_temp)
            sampling_temp="$2"
            shift 2
            ;;
        --model_name)
            model_name="$2"
            shift 2
            ;;
        --seed)
            seed="$2"
            shift 2
            ;;
        --limit)
            limit="$2"
            shift 2
            ;;
        --device)
            device="$2"
            shift 2
            ;;
        --chain)
            chain="$2"
            shift 2
            ;;
        --ids)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                ids+=("$1")
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

if [ -z "$structs_dir" ]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Convert structs_dir to absolute path if relative
if [[ "$structs_dir" != /* ]]; then
    structs_dir="$(cd "$structs_dir" && pwd)"
fi
# Convert save_path to absolute path if relative
if [[ -n "$save_path" && "$save_path" != /* ]]; then
    save_path="$(cd "$(dirname "$save_path")" && pwd)/$(basename "$save_path")"
fi

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$BASH_SOURCE")
cd "$SCRIPT_DIR/.."
. ./paths.sh # Load PROTEINMPNN_ENV (= ESMFOLD_ENV: needs both ProteinMPNN + ESMFold)

eval "$(conda shell.bash hook)"
conda activate "$PROTEINMPNN_ENV"
# Fix for Karolina compute nodes whose /lib64/libstdc++.so.6 lacks GLIBCXX_3.4.29
# (required by the env's pandas/numpy C extensions). Prepend the env's own libstdc++.
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
echo "Active conda environment: $(conda info --json | python -c "import sys, json; print(json.load(sys.stdin)['active_prefix_name'])")"
echo "Using python: $(which python)"

cd src/structure_metrics

args=("$structs_dir")
if [[ -n "$save_path" ]]; then
    args+=(--save_path "$save_path")
fi
if [[ -n "$num_seqs" ]]; then
    args+=(--num_seqs "$num_seqs")
fi
if [[ -n "$sampling_temp" ]]; then
    args+=(--sampling_temp "$sampling_temp")
fi
if [[ -n "$model_name" ]]; then
    args+=(--model_name "$model_name")
fi
if [[ -n "$seed" ]]; then
    args+=(--seed "$seed")
fi
if [[ -n "$limit" ]]; then
    args+=(--limit "$limit")
fi
if [[ -n "$device" ]]; then
    args+=(--device "$device")
fi
if [[ -n "$chain" ]]; then
    args+=(--chain "$chain")
fi
if [[ ${#ids[@]} -gt 0 ]]; then
    args+=(--ids "${ids[@]}")
fi

python run_self_consistency.py "${args[@]}"
