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

ENZYME_EXPLORER_PATH="/home2/soldat/documents/terpene_synthases/EnzymeExplorer"
ENZYME_EXPLORER_ENV="enzyme_explorer" # Enzyme Explorer conda environment name
ENZYME_EXPLORER_SEQUENCE_ONLY_PATH=$ENZYME_EXPLORER_PATH
ENZYME_EXPLORER_SEQUENCE_ONLY_ENV=$ENZYME_EXPLORER_ENV # Enzyme Explorer (sequence only) conda environment name

############################################################
# Broad homology search (Swiss-Prot + AlphaFold-Swiss-Prot)#
############################################################
# Both searches reuse the tps_eval env (DIAMOND + foldseek live there). They
# classify each hit TPS vs non-TPS by membership in the committed accession set.

# TPS accession set — COMMITTABLE default (lives in the repo). Override only if you
# regenerate it elsewhere. Generated via the UniProt REST query:
#   (reviewed:true) AND ((ec:4.2.3.*) OR (ec:5.5.1.*))
TPS_ACCESSIONS="$TPS_EVAL_ROOT/src/homology_search/tps_uniprot_accessions.txt"

# DIAMOND DB built from uniprot_sprot.fasta (diamond makedb). PER-INSTALL absolute
# path — built ON the cluster OUTSIDE the repo, never committed. The value below is
# a placeholder to set per-install (like SOLUPROT_PATH).
SWISSPROT_DIAMOND_DB="/home/soldat/documents/databases/swissprot_diamond/swissprot"

# foldseek AlphaFold/Swiss-Prot DB (foldseek databases "Alphafold/Swiss-Prot" ...).
# PER-INSTALL absolute path — downloaded ON the cluster OUTSIDE the repo, never
# committed. Placeholder to set per-install.
AFDB_SWISSPROT_DB="/home/soldat/documents/databases/afdb_swissprot/afdb_swissprot"
