using Pkg

cd(@__DIR__)
# TODO create julia project in tps_eval
# Pkg.activate(".") # TODO uncomment
Pkg.activate("../../../../terpene_generation/src") # TODO delete
#Pkg.instantiate()

# --- Script ---------------------------------------------------------
include("min_embedding_distance.jl")

num_args = length(ARGS)

if num_args == 1
    embeddings_path = ARGS[1]
    min_embedding_distance(embeddings_path; save=true)
elseif num_args == 2
    embeddings_path = ARGS[1]
    train_embeddings_path = ARGS[2]
    min_embedding_distance(embeddings_path; train_embeddings_path=train_embeddings_path, save=true)
else
    error("Invalid number of arguments. Expected 1 or 2, got $num_args.")
end
