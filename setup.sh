#!/bin/bash

conda create -n tps_eval -c conda-forge -c bioconda -c schrodinger python biopython pandas requests tqdm openbabel foldseek pymol-bundle -y

conda activate tps_eval
pip install torch
pip install fair-esm
pip install -e .
