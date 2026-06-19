#!/bin/bash
#
# compute_reference_stats.sh — reference-statistics pipeline for tps_eval.
#
# STANDALONE tool, deliberately NOT wired into scripts/run_eval_pipeline.py (the
# submit-all orchestrator). It runs the APPLICABLE metric tools on the MARTS-DB
# known-TPS reference set and aggregates each metric column into summary
# statistics ("natural TPS" bands), producing the committable reference-stats
# JSON (src/reference_stats/marts_db_metric_stats.json) that the rest of the
# project loads to compare generated designs against the natural distribution.
#
# It does NOT reimplement any metric — it submits the existing per-tool SLURM job
# scripts (via scripts/submit_job.sh) on the MARTS-DB fasta / structures, collects
# the resulting per-design CSVs into a reference output dir OUTSIDE the repo (like
# the DBs), and chains the aggregation step as a dependent job.
#
# ---------------------------------------------------------------------------
# WHICH METRICS GET A NATURAL BAND (and why)
# ---------------------------------------------------------------------------
# DEFAULT = run/band ALL intrinsic metrics (a "natural TPS distribution" is
# meaningful), NOT a curated subset:
#   sequence : motif_pair_distance, esm_pseudo_perplexity, motif_search,
#              soluprot, enzyme_explorer_sequence_only
#   structure: plddt, motif_structural_distance, active_site_geometry,
#              aggregation, domain_composition, proteinmpnn_score,
#              radius_of_gyration, aromatic_lining, diphosphate_sensor,
#              pocket_descriptors, ion_site_check
#   PAE (needs --pae_dir): global_confidence, interdomain_pae
#   holo-only (separate build): substrate_positioning (all-NaN on apo)
# The aggregator (aggregate_reference_stats.py) bands EVERY metric CSV it finds,
# so a newly-added tool's CSV is picked up with no edit here.
# EXCLUDED — inherently COMPARATIVE ("similarity to a SEPARATE reference set"), so a
# natural band is not meaningful (the reference set IS the natural set):
#   max_sequence_identity, local_sequence_search, min_embedding_distance,
#   structural_identity, domain_structural_identity, swissprot_search,
#   foldseek_swissprot_search, knn_label_transfer, sdr_divergence, substrate_class
# (self_consistency / scRMSD is a design-faithfulness metric, also excluded.)
#
# ---------------------------------------------------------------------------
# USAGE
# ---------------------------------------------------------------------------
#   scripts/compute_reference_stats.sh --cluster aurum \
#       --fasta_path data/train/TPS_sequences.fasta \
#       --ref_dir /home/soldat/documents/databases/marts_db_reference_stats \
#       [--structs_dir <dir_of_MARTS-DB_structures>] \
#       [--pae_dir <dir_of_saved_PAE_npz>] \
#       [--sequence_only] [--aggregate_only]
#
# * --fasta_path   MARTS-DB known-TPS FASTA (the marts_E* records).
# * --pae_dir      Dir of saved PAE npz (enables global_confidence + interdomain_pae).
# * --ref_dir      Reference OUTPUT dir on the cluster, OUTSIDE the repo. The
#                  metric CSVs are collected here; the aggregator reads them.
# * --structs_dir  Dir of MARTS-DB structures (.pdb/.cif or an af_output tree).
#                  If omitted, only the SEQUENCE metrics run; the structure
#                  metrics are skipped (documented TODO — see Step-1 report).
# * --sequence_only  Force-skip structure metrics even if --structs_dir is given.
# * --aggregate_only Skip submission; just aggregate CSVs already in --ref_dir
#                    (run locally, no SLURM). Handy after the metric jobs finish.
#
# The metric tools write next to their input (<input>_<tool>.csv). This script
# stages a COPY of fasta/structs into --ref_dir so the CSVs land there, keeping
# the repo and the gitignored data tree clean.

set -euo pipefail

CLUSTER=""
FASTA_PATH=""
STRUCTS_DIR=""
PAE_DIR=""
REF_DIR=""
SEQUENCE_ONLY="false"
AGGREGATE_ONLY="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cluster)       CLUSTER="$2"; shift 2 ;;
        --fasta_path)    FASTA_PATH="$2"; shift 2 ;;
        --structs_dir)   STRUCTS_DIR="$2"; shift 2 ;;
        --pae_dir)       PAE_DIR="$2"; shift 2 ;;
        --ref_dir)       REF_DIR="$2"; shift 2 ;;
        --sequence_only) SEQUENCE_ONLY="true"; shift ;;
        --aggregate_only) AGGREGATE_ONLY="true"; shift ;;
        -h|--help)
            grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -z "$REF_DIR" ]]; then
    echo "ERROR: --ref_dir is required (reference output dir, OUTSIDE the repo)."
    exit 1
fi
mkdir -p "$REF_DIR"
REF_DIR="$(cd "$REF_DIR" && pwd)"

OUTPUT_JSON="$REPO_ROOT/src/reference_stats/marts_db_metric_stats.json"

# ---------------------------------------------------------------------------
# Aggregate-only fast path (no SLURM, runs locally in the tps_eval env).
# ---------------------------------------------------------------------------
if [[ "$AGGREGATE_ONLY" == "true" ]]; then
    echo "[aggregate-only] Aggregating CSVs in $REF_DIR -> $OUTPUT_JSON"
    sh "$SCRIPT_DIR/run_aggregate_reference_stats.sh" \
        --input_dir "$REF_DIR" --output "$OUTPUT_JSON" --reference_name marts_db
    exit 0
fi

if [[ -z "$CLUSTER" || -z "$FASTA_PATH" ]]; then
    echo "ERROR: --cluster and --fasta_path are required (unless --aggregate_only)."
    exit 1
fi
FASTA_PATH="$(cd "$(dirname "$FASTA_PATH")" && pwd)/$(basename "$FASTA_PATH")"

# Stage a copy of the reference fasta into REF_DIR so the sequence-metric CSVs
# (named <fasta>_<tool>.csv) land in REF_DIR, not next to the source fasta.
REF_FASTA="$REF_DIR/$(basename "$FASTA_PATH")"
cp -f "$FASTA_PATH" "$REF_FASTA"
echo "[stage] reference fasta -> $REF_FASTA"

submit() { # submit <job_name> <job_args...> ; echoes the job id
    local job_name="$1"; shift
    local out
    out="$(sh "$SCRIPT_DIR/submit_job.sh" --cluster "$CLUSTER" \
            --job_name "$job_name" --job_args "$@")"
    echo "$out" >&2
    echo "$out" | awk '/[0-9]+$/{id=$NF} END{print id}'
}

declare -a METRIC_JOB_IDS=()

echo "=== SEQUENCE metrics (MARTS-DB reference) ==="
# All INTRINSIC sequence metrics (a "natural TPS" distribution is meaningful). The
# comparative ones (max_sequence_identity, local_sequence_search, min_embedding_distance,
# swissprot_search) are intentionally NOT run/banded -- see the header.
METRIC_JOB_IDS+=("$(submit motif_pair_distance --fasta_path "$REF_FASTA")")
METRIC_JOB_IDS+=("$(submit esm_pseudo_perplexity --fasta_path "$REF_FASTA")")
METRIC_JOB_IDS+=("$(submit motif_search --fasta_path "$REF_FASTA")")
METRIC_JOB_IDS+=("$(submit soluprot --fasta_path "$REF_FASTA")")
METRIC_JOB_IDS+=("$(submit enzyme_explorer_sequence_only --fasta_path "$REF_FASTA")")

# ---------------------------------------------------------------------------
# STRUCTURE metrics — only if a MARTS-DB structures dir was supplied.
# ---------------------------------------------------------------------------
if [[ "$SEQUENCE_ONLY" != "true" && -n "$STRUCTS_DIR" ]]; then
    STRUCTS_DIR="$(cd "$STRUCTS_DIR" && pwd)"
    echo "=== STRUCTURE metrics (MARTS-DB reference; structs=$STRUCTS_DIR) ==="
    # The structure tools write <structs_dir>_<tool>.csv (sibling of the dir).
    # We point --save_path into REF_DIR for each so all CSVs collect there.
    base="$(basename "$STRUCTS_DIR")"
    METRIC_JOB_IDS+=("$(submit plddt --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_plddt.csv")")
    METRIC_JOB_IDS+=("$(submit motif_structural_distance --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_motif_structural_distance.csv")")
    METRIC_JOB_IDS+=("$(submit active_site_geometry --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_active_site_geometry.csv")")
    METRIC_JOB_IDS+=("$(submit aggregation --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_aggregation.csv")")
    METRIC_JOB_IDS+=("$(submit domain_composition --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_domain_composition.csv")")
    METRIC_JOB_IDS+=("$(submit proteinmpnn_score --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_proteinmpnn_score.csv")")
    METRIC_JOB_IDS+=("$(submit radius_of_gyration --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_radius_of_gyration.csv")")
    METRIC_JOB_IDS+=("$(submit aromatic_lining --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_aromatic_lining.csv")")
    METRIC_JOB_IDS+=("$(submit diphosphate_sensor --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_diphosphate_sensor.csv")")
    METRIC_JOB_IDS+=("$(submit pocket_descriptors --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_pocket_descriptors.csv")")
    # ion_site_check carries signal only on holo folds (modelled ions); on apo
    # structures it reports n_ions_modelled=0 -- harmless to band. substrate_positioning
    # is holo-ONLY (all-NaN on apo), so it is left to the holo (AF3-cofold / Boltz2) build.
    METRIC_JOB_IDS+=("$(submit ion_site_check --structs_dir "$STRUCTS_DIR" \
        --save_path "$REF_DIR/${base}_ion_site_check.csv")")
    # PAE-derived fold confidence: only if a saved-PAE dir is supplied.
    if [[ -n "$PAE_DIR" ]]; then
        METRIC_JOB_IDS+=("$(submit global_confidence --structs_dir "$STRUCTS_DIR" \
            --pae_dir "$PAE_DIR" --save_path "$REF_DIR/${base}_global_confidence.csv")")
        METRIC_JOB_IDS+=("$(submit interdomain_pae --structs_dir "$STRUCTS_DIR" \
            --pae_dir "$PAE_DIR" --save_path "$REF_DIR/${base}_interdomain_pae.csv")")
    else
        echo "[note] global_confidence + interdomain_pae skipped (no --pae_dir)."
    fi
else
    echo "=== STRUCTURE metrics SKIPPED (no --structs_dir) ==="
    echo "    Sequence-metric stats only. To add structure stats, supply a dir"
    echo "    of MARTS-DB structures via --structs_dir (see Step-1 report:"
    echo "    download via src/alphafold/alphafold_struct_downloader.py keyed"
    echo "    by src/homology_search/tps_uniprot_accessions.txt)."
fi

# Build the afterok dependency list of all submitted metric jobs.
DEP=""
for jid in "${METRIC_JOB_IDS[@]}"; do
    [[ -n "$jid" ]] && DEP="${DEP:+$DEP:}$jid"
done
echo "Submitted metric jobs: $DEP"

# ---------------------------------------------------------------------------
# Chain the aggregation step after all metric jobs succeed.
# ---------------------------------------------------------------------------
if [[ -f "$SCRIPT_DIR/$CLUSTER/jobs/aggregate_reference_stats.sh" ]]; then
    echo "=== Aggregation (dependent on metric jobs) ==="
    sh "$SCRIPT_DIR/submit_job.sh" --cluster "$CLUSTER" \
        --job_name aggregate_reference_stats \
        --submit_args "--dependency=afterok:$DEP" \
        --job_args --input_dir "$REF_DIR" --output "$OUTPUT_JSON" \
            --reference_name marts_db
else
    echo "[note] No SLURM job script for aggregate_reference_stats on $CLUSTER."
    echo "       After the metric jobs finish, aggregate with:"
    echo "         scripts/compute_reference_stats.sh --ref_dir $REF_DIR --aggregate_only"
fi
