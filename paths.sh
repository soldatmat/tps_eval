#!/bin/bash

TPS_EVAL_ROOT="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"

############################################################
# Installation-depended paths & conda enviroment names     #
############################################################

TPS_EVAL_ENV="tps_eval" # tps_eval conda environment name

SOLUPROT_PATH="/home2/soldat/documents/soluprot"
SOLUPROT_ENV="soluprot" # SoluProt conda environment name

ENZYME_EXPLORER_PATH="/home2/soldat/documents/terpene_synthases/EnzymeExplorer"
ENZYME_EXPLORER_ENV="enzyme_explorer" # Enzyme Explorer conda environment name
ENZYME_EXPLORER_SEQUENCE_ONLY_PATH=$ENZYME_EXPLORER_PATH
ENZYME_EXPLORER_SEQUENCE_ONLY_ENV=$ENZYME_EXPLORER_ENV # Enzyme Explorer (sequence only) conda environment name
