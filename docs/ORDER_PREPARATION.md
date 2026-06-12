# Order preparation

Turn amino-acid designs into synthesis-ready DNA for **Golden Gate / MoClo** plasmid
construction: codon-optimize each protein, then add the fixed flanking overhangs.

This is a **standalone ordering utility**, *not* part of the evaluation / submit-all
pipeline — it is not an orchestrator Step, does not emit an `ID`-keyed metric CSV, and
needs no SLURM (it runs on a login node in seconds). Logic lives in
[`src/order_preparation/`](../src/order_preparation); the wrapper is
[`scripts/run_prepare_order.sh`](../scripts/run_prepare_order.sh).

## Pipeline (per design)

1. **Codon-optimize** the protein into a CDS for the target organism, *removing internal
   BsaI/BsmBI sites* — [`codon_optimization.py`](../src/order_preparation/codon_optimization.py)
2. **Add the fixed overhangs** around the CDS (the flanks are never modified) —
   [`overhangs.py`](../src/order_preparation/overhangs.py)
3. **Validate** the assembled construct — [`prepare_order.py`](../src/order_preparation/prepare_order.py)

## Usage

```sh
# Batch from a FASTA (default: yeast, Type 3 overhangs):
bash scripts/run_prepare_order.sh designs.fasta

# From a CSV, with a custom output prefix:
bash scripts/run_prepare_order.sh designs.csv -o my_order --organism yeast

# A single inline sequence (prints "id,sequence" to stdout):
bash scripts/run_prepare_order.sh --sequence MGRSY...PIPL --id design1
```

| Argument | Default | Meaning |
|----------|---------|---------|
| `input_path` | — | FASTA or CSV/TSV of protein designs (required unless `--sequence`) |
| `--sequence` | — | A single amino-acid sequence instead of an input file |
| `--id` | `design` | ID for `--sequence` mode |
| `-o`, `--output_prefix` | input without extension | Output path prefix |
| `--organism` | `yeast` | Target organism for codon usage (see below) |
| `--overhang_type` | `Type 3` | Golden Gate overhang type (see below) |
| `--seq-column` | auto-detect | Amino-acid column name (CSV input) |
| `--id-column` | auto-detect | ID column name (CSV input) |
| `--method` | `match_codon_usage` | DNAChisel CodonOptimize method (`match_codon_usage` \| `use_best_codon` \| `harmonize_rca`) |
| `--max_homopolymer` | `6` | Longest allowed single-nucleotide run (0 disables) |
| `--gc_min` / `--gc_max` | `0.30` / `0.65` | GC-window bounds (fraction) |
| `--gc_window` | `50` | GC sliding-window size in bp (0 disables) |
| `--seed` | `0` | RNG seed for reproducible sampling (`-1` = nondeterministic) |
| `--max_attempts` | `5` | Re-optimization tries to clear a Type IIS site before failing a design |

## Input

Auto-detected by extension:
- **FASTA** (`.fasta/.fa/.faa/.fas/.fna`) — the record ID is the construct ID.
- **CSV/TSV** (`.csv/.tsv`) — the ID and amino-acid columns are found by name
  (`id`/`name`/… and `sequence`/`aa`/`protein`/…), or set explicitly with
  `--id-column` / `--seq-column`. With no ID column, the row index is used.

## Output

- `<prefix>_order.csv` — one row per design: `id, status, protein, cds, ordered_sequence,
  length_nt, warnings`. Includes **every** design, even `FAILED` ones (for inspection).
- `<prefix>_order.txt` — `id,ordered_sequence` per line, matching the format of the
  previous order file (`all_candidates_dna_fixed.txt`) for direct submission. Contains
  **only `status == "ok"` designs** — any `FAILED` design is excluded.

`status` is `ok` or `FAILED` (a Type IIS site that could not be cleared — see below).
`warnings` is empty when all soft checks pass; any non-empty value is also printed to the
console (`[WARN]` / `[FAIL]`). The run **exits non-zero if any design FAILED**, so a broken
construct is never silently submitted.

## Codon optimization

Backend: **[DNAChisel](https://edinburgh-genome-foundry.github.io/DnaChisel/)** + codon
tables from `python_codon_tables` (both pinned in
[`requirements.txt`](../src/order_preparation/requirements.txt)). For each protein it
reverse-translates, enforces the translation, and optimizes synonymous codons for the
target organism subject to three **hard constraints** (codon usage is a soft objective
optimized within them):

1. **No internal BsaI (`GGTCTC`) / BsmBI (`CGTCTC`) sites** (both strands) — those Type IIS
   sites are introduced *only* by the overhangs; a copy inside the CDS would fragment the
   part during assembly.
2. **Homopolymer cap** (`--max_homopolymer`, default 6) — no single-nucleotide run longer
   than 6; long runs cause synthesis slippage/indels and can mimic cryptic yeast termination.
3. **GC window** (`--gc_min/--gc_max/--gc_window`, default 30–65 % over 50 bp) — keeps
   *local* GC in range, what synthesis vendors actually check.

Constraints 2–3 apply to **both** codon methods, so the output is synthesis-clean either
way. On the 32-design previous batch, defaults give mean GC 40 % and max homopolymer 6 —
matching the cad-sge web tool's profile while *also* guaranteeing the Type-IIS cleanliness
the web tool lacked.

- **Method** — `--method` default `match_codon_usage` samples codons to match the
  organism's natural distribution (like the web tool: natural GC + codon diversity).
  `use_best_codon` maximizes CAI but is lower-GC and repetitive (still synthesis-clean here
  thanks to the constraints, just less diverse). Sampling is seeded (`--seed`, default 0)
  so a design reproducibly yields the same sequence; pass `--seed -1` for fresh draws.
- **Graceful relaxation** — if the homopolymer/GC constraints make a particular protein
  infeasible (an amino-acid stretch with only AT-rich synonymous codons can fight a GC
  floor), they are loosened step by step (widen GC window → drop it → raise/drop the
  homopolymer cap) until a solution is found, and a note is added to that design's
  `warnings`. **Translation and enzyme-site avoidance are never relaxed.**
- **Organism** — default `yeast` → *Saccharomyces cerevisiae* (`s_cerevisiae_4932`). A
  small alias map (`yeast`, `s_cerevisiae`, `e_coli`, `human`, …) resolves friendly names;
  any `python_codon_tables` identifier/taxid also works directly, so adding organisms is a
  one-line change.
- **Stop codon** — appended as part of the optimized CDS, so the optimizer picks the
  organism-preferred stop (TAA in yeast, matching 23/32 of the previous order). The flanks
  are added *after* and never re-optimized.

> Why not the cad-sge.com web tool you may have used before? It is the open-source iDOG
> toolkit — bacterially oriented (RBS/terminator design, heavy ViennaRNA/TranstermHP deps),
> has no API, and does **not** remove Type IIS sites. DNAChisel installs locally, matches
> the web tool's yeast GC/codon profile, and additionally guarantees Type-IIS-clean CDSs
> and the homopolymer/GC limits above.

## Overhangs

Transcribed from `Plasmid_Generator.xlsx`. Every type shares a constant scaffold and
differs only in its fusion overhangs:

```
5' :  actcgacaac | CGTCTCa (BsmBI) | tcGGTCTCa (BsaI) | <left fusion>  … CDS …
3' :  <right fusion> | tGAGACC (BsaI rev) | tGAGACG (BsmBI rev) | gttgtggtgt
```

The default and only **validated** type is **`Type 3`** ("heterologous protein
expression"), whose flanks were cross-checked byte-exact against all 32 sequences of a
previous real order:

```
prefix = actcgacaacCGTCTCatcGGTCTCaT          (then ATG…)
suffix = ATCCtGAGACCtGAGACGgttgtggtgt         (preceded by the CDS stop codon)
```

All 16 sheet types are encoded for future use, but only `Type 3` is validated. (Sheet
note: *"Type 3 overhangs are only T on purpose. Don't touch them."*) The CDS-type rows
bake a `TAG` stop into the right-fusion cell; it is stripped at assembly because the stop
comes from codon optimization.

## Validation

`validate_construct` independently re-checks each assembled sequence (`EnforceTranslation`
already guarantees the protein, but this is a cheap safety net): starts with `ATG` after
the prefix, ends in a stop codon before the suffix, CDS length is a multiple of 3, the CDS
translates back to the input protein, **no BsaI/BsmBI site lies inside the CDS or across a
junction** (only the deliberate flank sites are allowed), **no homopolymer run exceeds the
cap**, and **the CDS GC stays within the window**.

Findings are split by severity:

- **Type IIS site (FATAL)** — a BsaI/BsmBI site inside the CDS or spanning a junction would
  be cut during Golden Gate, so it is *never* shipped. The optimizer's hard constraint keeps
  the CDS clean; a junction site (rare, depends on the sampled codons) triggers **automatic
  re-optimization with a fresh seed**, up to `--max_attempts` (default 5). If it still
  persists, the design is marked `status = FAILED`, **excluded from the order `.txt`**, kept
  in the CSV for inspection, and the run **exits non-zero**. *(The previous manual web
  workflow let 2 of 32 sequences keep internal BsaI sites — this prevents that class of
  defect from ever reaching an order.)*
- **Homopolymer / GC (SOFT)** — recorded in `warnings` but never fail a design. Because these
  checks mirror the optimization targets, a sequence whose constraints had to be relaxed
  flags here too, so the relaxation note and the residual reality both reach the user.

## Dependencies

```sh
pip install -r src/order_preparation/requirements.txt   # into the tps_eval conda env
```

`dnachisel` + `python_codon_tables`; `biopython` and `pandas` already ship with the
`tps_eval` env.
