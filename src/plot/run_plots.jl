using Pkg

cd(@__DIR__)
# TODO create julia project in tps_eval
Pkg.activate("../../../../terpene_generation/src") # TODO delete
#Pkg.instantiate()

# --- Script ---------------------------------------------------------
include("plot_comparison.jl")

num_args = length(ARGS)

# Parse required arguments
num_args < 3 && error("Invalid number of arguments. Expected at least 3: fasta_paths, data_names, data_colors.")
fasta_paths = Vector{String}(split(ARGS[1], ","))
data_names = Vector{String}(split(ARGS[2], ","))
data_colors = Vector{Symbol}(Symbol.(split(ARGS[3], ",")))

# Parse optional arguments
kwargs = Dict{Symbol,Any}()
i = 4
while i <= num_args
    if ARGS[i] == "--targets"
        if i + 1 > num_args
            error("Missing value for --targets.")
        end
        kwargs[:targets] = Vector{String}(split(ARGS[i+1], ","))
        global i += 2
    elseif ARGS[i] == "--save_dir"
        if i + 1 > num_args
            error("Missing value for --save_dir.")
        end
        kwargs[:save_dir] = ARGS[i+1]
        global i += 2
    else
        error("Unknown argument: $(ARGS[i])")
    end
end

plot_comparison(
    fasta_paths,
    data_names,
    data_colors;
    kwargs...
)
