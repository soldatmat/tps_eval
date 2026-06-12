#!/bin/bash

git submodule update --init --recursive

conda create -n tps_eval -c conda-forge -c bioconda python biopython pandas matplotlib scipy scikit-learn requests tqdm openbabel foldseek diamond mmseqs2 pymol-open-source -y

conda activate tps_eval
pip install torch
pip install fair-esm
pip install -e .
