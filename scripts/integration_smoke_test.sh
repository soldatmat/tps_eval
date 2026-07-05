#!/usr/bin/env bash
# Light end-to-end integration smoke test for tps_eval.
#
# Unlike the per-module unit tests (src/**/test_*.py, which use synthetic data and
# mock external tools), this drives the REAL run_<tool>.py argv entrypoints on tiny
# real inputs and verifies each writes a well-formed, ID-keyed CSV. It covers only
# the tools that need NO external binary / model / GPU / reference DB, so it runs in
# seconds inside an activated tps_eval env on any cluster (or a laptop):
#
#     conda activate tps_eval        # numpy, pandas, biopython, dnachisel
#     bash scripts/integration_smoke_test.sh
#
# Exits non-zero if any entrypoint fails or produces a malformed/empty output.
set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$REPO_ROOT/src"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

pass=0; fail=0; failed=()

# --- build tiny real inputs ---------------------------------------------------
# A FASTA with one class-I TPS-like sequence carrying a DDXXD + NSE motif, and one
# plain sequence (so both a hit row and a NaN row are exercised).
cat > "$WORK/designs.fasta" <<'FASTA'
>design_1
MSTLPISKVDDIYDVYGDKSEEILAFTRAFDRWDVNSEKTLPEYMKMAFASLYNFVNEHAY
>design_2
MKAILVGGAGYIGSHTVVELLEAGHEVVVVDNLSNGHREAVPKGVPFYEGDIRDRALLDRV
FASTA

# Two minimal structures (Cα-only), B-factor field carries a pLDDT-like value so
# run_plddt reads it; coordinates give run_radius_of_gyration a finite Rg.
python - "$WORK" <<'PY'
import os, sys
work = sys.argv[1]
structs = os.path.join(work, "structs")
os.makedirs(structs, exist_ok=True)
def write_pdb(path, coords, bfac):
    with open(path, "w") as fh:
        for i, (x, y, z) in enumerate(coords, start=1):
            fh.write(
                f"ATOM  {i:>5d}  CA  ALA A{i:>4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{bfac:6.2f}           C\n"
            )
        fh.write("END\n")
# 4 Cα spread along a line; distinct B-factors.
write_pdb(os.path.join(structs, "design_1.pdb"),
          [(0,0,0),(3.8,0,0),(7.6,0,0),(11.4,0,0)], 88.0)
write_pdb(os.path.join(structs, "design_2.pdb"),
          [(0,0,0),(0,3.8,0),(0,7.6,0),(0,11.4,0)], 72.5)
PY

# A single short, cleanly-orderable design for the batch order-preparation path
# (--sequence mode prints to stdout by design; only batch mode writes a CSV).
cat > "$WORK/order_input.fasta" <<'FASTA'
>design_1
MSTLPISKVDDIYDVYGDKSEEIL
FASTA

FASTA="$WORK/designs.fasta"
STRUCTS="$WORK/structs"

# --- helper: run an entrypoint (from its src subdir) + validate its CSV -------
check_csv() {  # <label> <csv_path> [min_rows]
  local label="$1" csv="$2" min="${3:-1}"
  python - "$csv" "$min" <<'PY'
import sys, pandas as pd
csv, need = sys.argv[1], int(sys.argv[2])
df = pd.read_csv(csv)
# Metric tools key by "ID"; the standalone order sheet keys by lowercase "id".
cols = {c.lower() for c in df.columns}
assert "id" in cols, f"missing id/ID column in {csv}: {list(df.columns)}"
assert len(df) >= need, f"{csv}: {len(df)} rows < {need}"
PY
}

run_case() {  # <label> <subdir> <expected_csv> -- <python argv...>
  local label="$1" subdir="$2" csv="$3"; shift 3; shift  # drop the literal --
  local out
  out=$( cd "$SRC/$subdir" && python "$@" 2>&1 )
  local rc=$?
  if [ $rc -eq 0 ] && [ -f "$csv" ] && check_csv "$label" "$csv" 2>/dev/null; then
    echo "PASS  $label"; pass=$((pass+1))
  else
    echo "FAIL  $label"; echo "$out" | tail -12 | sed 's/^/      | /'
    [ -f "$csv" ] || echo "      | (no output CSV at $csv)"
    fail=$((fail+1)); failed+=("$label")
  fi
}

echo "== tps_eval integration smoke test =="

# Sequence branch (no external deps): motif search + motif-pair distance.
run_case "motif_search"        sequence_metrics "${FASTA%.fasta}_motifs.csv" \
    -- run_motif_search.py "$FASTA" DDXXD "NSE/DTE"
run_case "motif_pair_distance" sequence_metrics "${FASTA%.fasta}_motif_pair_distance.csv" \
    -- run_motif_pair_distance.py "$FASTA"

# Structure branch (biopython only): radius of gyration + pLDDT.
run_case "radius_of_gyration"  structure_metrics "${STRUCTS}_radius_of_gyration.csv" \
    -- run_radius_of_gyration.py "$STRUCTS"
run_case "plddt"               structure_metrics "${STRUCTS}_plddt.csv" \
    -- run_plddt.py "$STRUCTS"

# Standalone order-preparation, batch mode (real dnachisel codon optimization →
# Golden Gate overhangs → <prefix>_order.csv).
run_case "prepare_order"       order_preparation "$WORK/order_order.csv" \
    -- run_prepare_order.py "$WORK/order_input.fasta" -o "$WORK/order"

echo "======================================"
echo "PASS=$pass  FAIL=$fail"
if [ "$fail" -gt 0 ]; then printf 'FAILED: %s\n' "${failed[@]}"; exit 1; fi
