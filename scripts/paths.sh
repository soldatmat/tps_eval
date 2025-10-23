#!/bin/bash

TPS_EVAL_ROOT="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"

############################################################
# Installation-depended paths & conda enviroment names     #
############################################################

TPS_EVAL_ENV="tps_eval" # tps_eval conda environment name

SOLUPROT_PATH="/home2/soldat/documents/soluprot"
SOLUPROT_ENV="soluprot" # SoluProt conda environment name

ENZYME_EXPLORER_SEQUENCE_ONLY_PATH="/home2/soldat/documents/TerpeneMiner"
ENZYME_EXPLORER_SEQUENCE_ONLY_ENV="terpene_miner" # Enzyme Explorer (sequence only) conda environment name
ENZYME_EXPLORER_PATH="/home2/soldat/documents/TerpeneMiner_easy"
ENZYME_EXPLORER_ENV="terpene_miner" # Enzyme Explorer conda environment name
