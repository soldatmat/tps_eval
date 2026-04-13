using CSV
using DataFrames

include("../data/sequences.jl")

convert_to_regex(motif::String)::Regex = Regex(motif)
convert_motifs(motifs::AbstractVector{String})::AbstractVector{Regex} = map(convert_to_regex, motifs)

function load_df(fasta_path::String)
    sequence_identifiers, sequences = separate_identifiers(load_fasta_sequences(fasta_path; load_identifiers=true))
    return DataFrame(ID=sequence_identifiers, sequence=sequences)
end

function find_motifs(
    df::DataFrame,
    motifs::AbstractVector{Regex},
)
    for motif in motifs
        motif_present = map(row -> occursin(motif, row.sequence), eachrow(df))
        column_name = replace(string(motif), r"^r\"|\"$" => "")
        df[!, column_name] = motif_present
    end
end

function get_save_path(data_path::String)
    extension = split(data_path, ".")[end]
    save_path = data_path[1:end-length(extension)-1] * "_motifs.csv"
    return save_path
end

########## Main ##########
function motif_search(
    fasta_path::String,
    motifs::AbstractVector{String};
    save::Bool=true,
)
    motifs = convert_motifs(motifs)
    motif_search(fasta_path, motifs; save=save)
end

function motif_search(
    fasta_path::String,
    motifs::AbstractVector{Regex};
    save::Bool=true,
)
    df = load_df(fasta_path)
    find_motifs(df, motifs)

    if save
        save_path = get_save_path(fasta_path)
        CSV.write(save_path, df)
    end

    return df
end
