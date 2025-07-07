using StatsBase
using CairoMakie

include("constants.jl")
include("../data/load_results.jl")


"""
    Expects all data to be saved in the same folder as the fasta file for each dataset. See `load_results` for details.
"""
function density_comparison(
    fasta_paths::Vector{String},
    data_names::Vector{String},
    data_colors::Vector{Symbol},
    target::String;
    direction::Symbol=:vertical, # :horizontal or :vertical
    save_dir::Union{String, Nothing}=nothing,
)
    ########## Get plot constans ##########
    load_list = LOAD[target]
    min_val = MIN_VAL[target]
    max_val = MAX_VAL[target]
    ticks = TICKS[target]
    threshold = THRESHOLD[target]
    offset = OFFSET[target]

    
    ########## Load data ##########
    all_dfs = map(fasta_path -> load_results(fasta_path; load=load_list), fasta_paths)
    all_data = map(df -> convert(Vector{Float64}, df[!, target]), all_dfs)


    ########## Density plot ##########

    if direction == :horizontal
        ########## Horizontal ##########
        fig = Figure()
        ax = Axis(
            fig[1, 1],
            title=target * " density",
            xticks=ticks,
            yticks=((1:length(data_names)) * offset, data_names),
        )

        for i in length(all_data):-1:1
            density!(
                ax,
                all_data[i],
                offset=offset * i,
                #boundary=(0.0, maximum(vcat(all_data...))),

                color=data_colors[i],
                #color=:x,
                #colormap=:thermal,
                #colorrange=(-5, 5),

                strokewidth=1,
                strokecolor=:black,
            )
        end
        xlims!(ax, low=min_val, high=max_val)

        isnothing(threshold) || vlines!(ax, threshold, color=:red, linestyle=:dash, linewidth=2)

        display(fig)

    elseif direction == :vertical
        ########## Vertical ##########
        fig = Figure()
        ax = Axis(
            fig[1, 1],
            title=target * " density",
            yticks=ticks,
            xticks=((1:length(data_names)) * offset, data_names),
        )

        for i in length(all_data):-1:1
            density!(
                ax,
                all_data[i],
                offset=offset * i,
                direction=:y,
                #boundary=(0.0, maximum(vcat(all_data...))),

                color=data_colors[i],
                #color=:x,
                #colormap=:thermal,
                #colorrange=(-5, 5),

                strokewidth=1,
                strokecolor=:black,
            )
        end
        ylims!(ax, low=min_val, high=max_val)

        isnothing(threshold) || hlines!(ax, threshold, color=:red, linestyle=:dash, linewidth=2)

        display(fig)
    end


    ########## Save plot ##########
    if !isnothing(save_dir)
        isdir(save_dir) || mkpath(save_dir)
        save(joinpath(save_dir, target * "_density_" * String(direction) * ".png"), fig)
    end 
end
