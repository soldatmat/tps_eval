#!/bin/bash

git submodule update --init --recursive

conda create -n tps_eval -c conda-forge -c bioconda python biopython pandas matplotlib scipy scikit-learn requests tqdm openbabel foldseek diamond mmseqs2 pymol-open-source -y

conda activate tps_eval
pip install torch
pip install fair-esm
pip install -e .
# Extra deps for the standalone order-preparation tool (dnachisel + codon tables).
pip install -r src/order_preparation/requirements.txt
