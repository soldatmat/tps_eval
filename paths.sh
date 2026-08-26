#!/bin/bash

TPS_EVAL_ROOT="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"

############################################################
# Installation-depended paths & conda enviroment names     #
############################################################

TPS_EVAL_ENV="tps_eval" # tps_eval conda environment name

ESMFOLD_ENV="esmfold" # ESMFold (structure prediction) conda environment name

# ProteinMPNN (sequence design / scoring; vendored at vendor/ProteinMPNN) reuses
# the ESMFold env by default — it is lightweight (small model, torch + numpy) and
# the self-consistency tool needs both ProteinMPNN and ESMFold in one env. Point
# at a dedicated env only if ProteinMPNN conflicts with the esmfold env.
PROTEINMPNN_ENV="$ESMFOLD_ENV" # ProteinMPNN conda environment name

AGGRESCAN3D_ENV="aggrescan3d" # Aggrescan3D (A3D, structure-based aggregation propensity) conda environment name

# Active-site pocket descriptors (fpocket geometric + P2Rank ML ligandability).
# Dedicated env created via conda-forge:
#   conda create -n pocket -c conda-forge fpocket openjdk=21 python=3.11 pandas numpy biopython
# (openjdk is required by P2Rank's `prank` launcher; the env puts `java` on PATH.)
POCKET_ENV="pocket" # fpocket + P2Rank conda environment name

# P2Rank prebuilt distribution dir (must contain the `prank` launcher). PER-INSTALL
# absolute path — the official release tarball unpacked ON the cluster OUTSIDE the
# repo, never committed (like SOLUPROT_PATH). Download from:
#   https://github.com/rdk/p2rank/releases  (p2rank_<ver>.tar.gz)
# Leave empty / unset to skip the P2Rank cross-check (fpocket still runs).
P2RANK_PATH="/home/soldat/documents/tools/p2rank_2.5.1"

SOLUPROT_PATH="/home2/soldat/documents/soluprot"
SOLUPROT_ENV="soluprot" # SoluProt conda environment name

# CataPro (enzyme kinetics; vendored at vendor/CataPro) and TmProt (melting temperature;
# vendored at vendor/TmProt) each run in their own conda env, created by
# scripts/setup/setup_catapro.sh / scripts/setup/setup_tmprot.sh. Env NAMES only — the code is vendored
# in-repo, so (unlike SoluProt) there is no external _PATH to set.
CATAPRO_ENV="catapro" # CataPro conda environment name
TMPROT_ENV="tmprot"   # TmProt conda environment name

# NOTE: $HOME is /home/soldat since the Aurum3 rebuild -- /home2 no longer exists.
ENZYME_EXPLORER_PATH="$HOME/documents/terpene_synthases/EnzymeExplorer"
ENZYME_EXPLORER_ENV="enzyme_explorer_prod" # Enzyme Explorer conda environment name
ENZYME_EXPLORER_SEQUENCE_ONLY_PATH=$ENZYME_EXPLORER_PATH
ENZYME_EXPLORER_SEQUENCE_ONLY_ENV=$ENZYME_EXPLORER_ENV # Enzyme Explorer (sequence only) conda environment name

# ⚠ EnzymeExplorer WITH STRUCTURES (`ee_struct`) uses a DIFFERENT checkout on purpose.
# The `revision` branch above HARDCODES `prefilter_pdbs_by_foldseek=True` in
# enzymeexplorer/src/prediction/domains.py and does not expose it on the CLI. That
# prefilter skips (query x template) USalign pairs with no plausible foldseek alignment --
# its own docstring admits a "small recall loss" -- and for a GATE a recall loss means
# falsely rejecting real designs. Measured: in the 2026-07-13 run only 1 of 56 structures
# passed the prefilter, so domain detection never actually ran on the other 55 and the
# resulting "0/56 domain hits" could not be distinguished from a prefilter artifact.
# The `main` branch defaults it to False AND exposes `--prefilter-pdbs-by-foldseek` as an
# opt-in, so pointing ee_struct at that checkout fixes this with no upstream patch.
# (A/B on matched knobs: revision and main produce bit-identical detection output --
# see projects/MARTS_domain_detections/README.md, "Vocabulary resolved".)
ENZYME_EXPLORER_STRUCT_PATH="$HOME/documents/terpene_synthases/EnzymeExplorer_main"
ENZYME_EXPLORER_STRUCT_ENV="enzyme_explorer_main"

############################################################
# Broad homology search (Swiss-Prot + AlphaFold-Swiss-Prot)#
############################################################
# Both searches reuse the tps_eval env (DIAMOND + foldseek live there). They
# classify each hit TPS vs non-TPS by membership in the committed accession set.

# TPS accession set — COMMITTABLE default (lives in the repo). Override only if you
# regenerate it elsewhere. Generated via the UniProt REST query:
#   (reviewed:true) AND ((ec:4.2.3.*) OR (ec:5.5.1.*))
TPS_ACCESSIONS="$TPS_EVAL_ROOT/src/tps_eval/homology_search/tps_uniprot_accessions.txt"

# DIAMOND DB built from uniprot_sprot.fasta (diamond makedb). PER-INSTALL absolute
# path — built ON the cluster OUTSIDE the repo, never committed. The value below is
# a placeholder to set per-install (like SOLUPROT_PATH).
SWISSPROT_DIAMOND_DB="/home/soldat/documents/databases/swissprot_diamond/swissprot"

# foldseek AlphaFold/Swiss-Prot DB (foldseek databases "Alphafold/Swiss-Prot" ...).
# PER-INSTALL absolute path — downloaded ON the cluster OUTSIDE the repo, never
# committed. Placeholder to set per-install.
AFDB_SWISSPROT_DB="/home/soldat/documents/databases/afdb_swissprot/afdb_swissprot"
