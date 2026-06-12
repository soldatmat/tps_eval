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
| `--method` | `use_best_codon` | DNAChisel CodonOptimize method |

## Input

Auto-detected by extension:
- **FASTA** (`.fasta/.fa/.faa/.fas/.fna`) — the record ID is the construct ID.
- **CSV/TSV** (`.csv/.tsv`) — the ID and amino-acid columns are found by name
  (`id`/`name`/… and `sequence`/`aa`/`protein`/…), or set explicitly with
  `--id-column` / `--seq-column`. With no ID column, the row index is used.

## Output

- `<prefix>_order.csv` — one row per design: `id, protein, cds, ordered_sequence,
  length_nt, warnings`.
- `<prefix>_order.txt` — `id,ordered_sequence` per line, matching the format of the
  previous order file (`all_candidates_dna_fixed.txt`) for direct submission.

`warnings` is empty when all validation checks pass; any non-empty value is also printed
to the console during the run.

## Codon optimization

Backend: **[DNAChisel](https://edinburgh-genome-foundry.github.io/DnaChisel/)** + codon
tables from `python_codon_tables` (both pinned in
[`requirements.txt`](../src/order_preparation/requirements.txt)). For each protein it
reverse-translates, enforces the translation, optimizes synonymous codons for the target
organism, and **avoids the BsaI (`GGTCTC`) and BsmBI (`CGTCTC`) recognition sites on both
strands** — those Type IIS sites are introduced *only* by the overhangs, and any copy
inside the CDS would fragment the part during assembly.

- **Organism** — default `yeast` → *Saccharomyces cerevisiae* (`s_cerevisiae_4932`). A
  small alias map (`yeast`, `s_cerevisiae`, `e_coli`, `human`, …) resolves friendly names;
  any `python_codon_tables` identifier/taxid also works directly, so adding organisms is a
  one-line change.
- **Stop codon** — appended as part of the optimized CDS, so the optimizer picks the
  organism-preferred stop (TAA in yeast, matching 23/32 of the previous order). The flanks
  are added *after* and never re-optimized.

> Why not the cad-sge.com web tool you may have used before? It is the open-source iDOG
> toolkit — bacterially oriented (RBS/terminator design, heavy ViennaRNA/TranstermHP deps),
> has no API, and does **not** remove Type IIS sites. DNAChisel installs locally, tracks the
> web tool's yeast output closely, and additionally guarantees Type-IIS-clean CDSs.

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
translates back to the input protein, and **no BsaI/BsmBI site lies inside the CDS or across
a junction** (only the deliberate flank sites are allowed). Failures populate the
`warnings` column. *(The previous manual web workflow let 2 of 32 sequences keep internal
BsaI sites — this step catches that class of defect.)*

## Dependencies

```sh
pip install -r src/order_preparation/requirements.txt   # into the tps_eval conda env
```

`dnachisel` + `python_codon_tables`; `biopython` and `pandas` already ship with the
`tps_eval` env.
