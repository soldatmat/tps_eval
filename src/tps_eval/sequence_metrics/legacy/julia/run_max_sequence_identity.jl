using Pkg

cd(@__DIR__)
# TODO create julia project in tps_eval
# Pkg.activate(".") # TODO uncomment
Pkg.activate("../../../../terpene_generation/src") # TODO delete
#Pkg.instantiate()

# --- Script ---------------------------------------------------------
include("max_sequence_identity.jl")

num_args = length(ARGS)

if num_args == 1
    fasta_path = ARGS[1]
    max_sequence_identity(fasta_path)
elseif num_args == 2
    fasta_path = ARGS[1]
    train_path = ARGS[2]
    max_sequence_identity(fasta_path; train_path=train_path)
else
    error("Invalid number of arguments. Expected 1 or 2, got $num_args.")
end
