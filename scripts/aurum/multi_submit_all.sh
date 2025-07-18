#!/bin/bash

models=(
    TPS_dplm_150m_stage3_run_7_best
    TPS_dplm_150m_stage3_run_8_step20000
    TPS_dplm_650m_stage3_run_1_best
    TPS_dplm_650m_stage3_run_3_step20000
)
temperatures=(0.0 1.0 4.0 8.0)

for model in "${models[@]}"; do
    for temp in "${temperatures[@]}"; do
        fasta_path="/home2/soldat/documents/terpene_synthases/output/dplm/${model}/sl1000_t${temp}/all_sequences.fasta"
        train_path="/home2/soldat/documents/terpene_synthases/data/MARTS-DB/2025-04-24/TPS_sequences.fasta"
        sh submit_all.sh --fasta_path "$fasta_path" --train_path "$train_path"
    done
done
