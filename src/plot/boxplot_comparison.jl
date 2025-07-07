using StatsBase
using CairoMakie

include("constants.jl")
include("../data/load_results.jl")


"""
    Expects all data to be saved in the same folder as the fasta file for each dataset. See `load_results` for details.

    Example:
    `boxplot_comparison(["/Users/soldatmat/Documents/terpene_synthases/tps_eval/data/generated/generated_sequences.fasta"], ["generated"], [:goldenrod1], "isTPS")`
"""
function boxplot_comparison(
    fasta_paths::Vector{String},
    data_names::Vector{String},
    data_colors::Vector{Symbol},
    target::String;
    save_dir::Union{String, Nothing}=nothing,
)
    ########## Get plot constans ##########
    load_list = LOAD[target]
    min_val = MIN_VAL[target]
    max_val = MAX_VAL[target]
    ticks = TICKS[target]
    threshold = THRESHOLD[target]

    
    ########## Load data ##########
    all_dfs = map(fasta_path -> load_results(fasta_path; load=load_list), fasta_paths)
    all_data = map(df -> convert(Vector{Float64}, df[!, target]), all_dfs)


    ########## Boxplot ##########
    categories = Vector{Int}(vcat([i * ones(length(all_data[i])) for i in eachindex(all_data)]...))
    data = vcat(all_data...)

    fig = Figure()
    ax = Axis(fig[1, 1];
        title=target * " boxplot",
        #xlabel="categories",
        xticks=(1:length(data_names), data_names),
        ygridvisible=true,
        ylabel=target,
        yticks=ticks,
        ygridcolor=:gray,
        ygridstyle=:dash,
    )
    ylims!(ax, low=min_val, high=max_val)

    boxplot!(ax, categories, data,
        color=map(c -> data_colors[c], categories),
        orientation=:vertical, # is default
    )

    isnothing(threshold) || hlines!(ax, threshold, color=:red, linestyle=:dash, linewidth=2)

    display(fig)


    ########## Save plot ##########
    if !isnothing(save_dir)
        isdir(save_dir) || mkpath(save_dir)
        save(joinpath(save_dir, target * "_boxplot.png"), fig)
    end 
end
