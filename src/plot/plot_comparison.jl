include("boxplot_comparison.jl")
include("density_comparison.jl")


function plot_comparison(
    fasta_paths::Vector{String},
    data_names::Vector{String},
    data_colors::Vector{Symbol},
    targets::Vector{String};
    save_dir::Union{String,Nothing}=nothing,
)
    for target in targets
        try
            println("Generating boxplot for target: $target...")
            boxplot_comparison(
            fasta_paths,
            data_names,
            data_colors,
            target;
            save_dir=save_dir,
            )

            println("Generating density plot for target: $target...")
            density_comparison(
            fasta_paths,
            data_names,
            data_colors,
            target;
            direction=:vertical,
            save_dir=save_dir,
            )
        catch e
            println("Error while plotting for target $target: $e")
            println("Stacktrace:")
            println(stacktrace(catch_backtrace()))
        end
    end
end

function main()
    FASTA_PATHS = [
        "/Users/soldatmat/Documents/terpene_synthases/tps_eval/data/generated/generated_sequences.fasta",
    ]

    DATA_NAMES = [
        "generated",
    ]

    DATA_COLORS = [
        :goldenrod1,
    ]

    TARGETS = [
        "sequence_identity",
        "sequence_identity_self",
        "sequence_similarity",
        "sequence_similarity_self",
        "min_embedding_distance",
        "min_embedding_distance_self",
        "isTPS",
        "isTPS_seq",
        "soluble",
    ]

    SAVE_DIR = "/Users/soldatmat/Documents/terpene_synthases/tps_eval/data/generated/plots/"

    plot_comparison(
        FASTA_PATHS,
        DATA_NAMES,
        DATA_COLORS,
        TARGETS;
        save_dir=SAVE_DIR,
    )
end
