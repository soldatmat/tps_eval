#!/bin/bash

USAGE="--fasta_paths <fasta1.fa> [<fasta2.fa> ...] --data_names <name1> [<name2> ...] --data_colors <color1> [<color2> ...] --targets <target1> [<target2> ...] [--save_dir <save_dir>]"

Help()
{
    # Display Help
    echo "Usage: $0 $USAGE"
    echo
    echo "Arguments:"
    echo "  --fasta_paths   One or more paths to FASTA files (required, space-separated list)"
    echo "  --data_names    One or more names corresponding to each FASTA file (required, space-separated list)"
    echo "  --data_colors   One or more colors corresponding to each dataset (required, space-separated list)"
    echo "  --targets       One or more targets corresponding to each dataset (required, space-separated list)"
    echo "  --save_dir      Directory to save output plots (optional)"
    echo "  -h, --help      Show this help message and exit"
    echo
    echo "Note: The number of elements for --fasta_paths, --data_names, --data_colors, and --targets must match."
}

# Parse long options manually
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_paths)
            shift
            fasta_paths=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                fasta_paths+=("$1")
                shift
            done
            ;;
        --data_names)
            shift
            data_names=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                data_names+=("$1")
                shift
            done
            ;;
        --data_colors)
            shift
            data_colors=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                data_colors+=("$1")
                shift
            done
            ;;
        --targets)
            shift
            targets=()
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                targets+=("$1")
                shift
            done
            ;;
        --save_dir)
            save_dir="$2"
            shift
            shift
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


if [[ -z "${fasta_paths[*]}" || -z "${data_names[*]}" || -z "${data_colors[*]}" || -z "${targets[*]}" ]]; then
    echo "Usage: $0 $USAGE"
    exit 1
fi

# Check that all four mandatory arguments have the same number of elements
if [[ ${#fasta_paths[@]} -ne ${#data_names[@]} || ${#fasta_paths[@]} -ne ${#data_colors[@]} ]]; then
    echo "Error: --fasta_paths, --data_names, and --data_colors must have the same number of elements."
    exit 1
fi

# Convert each fasta_path in fasta_paths to absolute path if it's relative
for i in "${!fasta_paths[@]}"; do
    if [[ "${fasta_paths[$i]}" != /* ]]; then
        fasta_paths[$i]="$(cd "$(dirname "${fasta_paths[$i]}")" && pwd)/$(basename "${fasta_paths[$i]}")"
    fi
done

############################################################
# Main                                                     #
############################################################
SCRIPT_DIR=$(dirname "$0")
cd "$SCRIPT_DIR/../src/plot"

if [[ -n "$save_dir" ]] && [[ "$save_dir" != "" ]]; then
    julia run_plots.jl "$fasta_paths" "$data_names" "$data_colors" "$targets" "$save_dir"
else
    julia run_plots.jl "$fasta_paths" "$data_names" "$data_colors" "$targets"
fi
