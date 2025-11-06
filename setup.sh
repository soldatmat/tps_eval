#!/bin/bash

conda create -n tps_eval -c conda-forge python biopython pandas requests tqdm -y

conda activate tps_eval
pip install torch
pip install fair-esm
