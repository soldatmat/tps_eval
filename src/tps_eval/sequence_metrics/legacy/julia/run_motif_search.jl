using Pkg

cd(@__DIR__)
# TODO create julia project in tps_eval
# Pkg.activate(".") # TODO uncomment
Pkg.activate("../../../../terpene_generation/src") # TODO delete
#Pkg.instantiate()

# --- Script ---------------------------------------------------------
include("motif_search.jl")

fasta_path = ARGS[1] # Path to the FASTA file with sequences
motifs = ARGS[2:end] # Motifs to search for, passed as command line arguments
motif_search(fasta_path, motifs; save=true)
