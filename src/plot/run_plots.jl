include("boxplot_comparison.jl")


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

for target in TARGETS
    println("Generating boxplot for target: $target...")
    boxplot_comparison(
        FASTA_PATHS,
        DATA_NAMES,
        DATA_COLORS,
        target;
        save_dir=SAVE_DIR,
    )
end
