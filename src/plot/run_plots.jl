using Pkg

cd(@__DIR__)
# TODO create julia project in tps_eval
Pkg.activate("../../../../terpene_generation/src") # TODO delete
#Pkg.instantiate()

# --- Script ---------------------------------------------------------
include("plot_comparison.jl")

num_args = length(ARGS)

if num_args == 4 || num_args == 5
    fasta_paths = split(ARGS[1], ",")
    data_names = split(ARGS[2], ",")
    data_colors = Symbol.(split(ARGS[3], ","))
    targets = split(ARGS[4], ",")
    save_dir = num_args == 5 ? ARGS[5] : nothing
else
    error("Invalid number of arguments. Expected 4 or 5, got $num_args.")
end

plot_comparison(
    fasta_paths,
    data_names,
    data_colors,
    targets;
    save_dir=save_dir
)
