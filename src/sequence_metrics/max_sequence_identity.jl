using CSV
using DataFrames
using BioSequences
using BioAlignments

include("../data/sequences.jl")
include("count_positives.jl")

SUBSTITUTION_MATRIX = BLOSUM62


function evaluate_max_sequence_identity(
    train,
    generated,
    generated_identifiers,
    train_identifiers;
    save_path::Union{String,Nothing}=nothing,
    return_second_max::Bool=false
)
    max_sequence_identity, max_sequence_similarity, max_sequence_identity_index, max_sequence_similarity_index = get_max_sequence_identity(train, generated; self_comparison=return_second_max)

    max_sequence_identity_hit = train_identifiers[max_sequence_identity_index]
    max_sequence_similarity_hit = train_identifiers[max_sequence_similarity_index]

    # Save calculated novelty into CSV
    if !isnothing(save_path)
        df = DataFrame(
            :ID => generated_identifiers,
            :sequence_identity => max_sequence_identity,
            :sequence_identity_hit => max_sequence_identity_hit,
            :sequence_similarity => max_sequence_similarity,
            :sequence_similarity_hit => max_sequence_similarity_hit,
        )
        CSV.write(save_path, df)
    end

    return max_sequence_identity, max_sequence_similarity, max_sequence_identity_hit, max_sequence_similarity_hit
end

"""
Computes sequence identity between each pair of sequences from `train` and `generated`
Returns max sequence identity for each `generated` sequence with any of the `train` sequences.

Sequence identity implementation:
- constructs pairwise alignment for each pair of sequences
- identity = #Matches / length(alignment)
"""
get_max_sequence_identity(sequences) = get_max_sequence_identity(sequences, sequences; self_comparison=true)

function get_max_sequence_identity(train, generated; self_comparison::Bool=false)
    # Convert sequences to ::LongAA
    train = LongAA.(train)
    generated = LongAA.(generated)
    get_max_sequence_identity(train, generated; self_comparison)
end

function get_max_sequence_identity(
    train::AbstractVector{LongAA},
    generated::AbstractVector{LongAA};
    self_comparison::Bool=false,
)
    # Construct pairwise alignments of each generated sequence to each train sequence
    score_model = AffineGapScoreModel(SUBSTITUTION_MATRIX, gap_open=-11, gap_extend=-1)
    align_sequences(seq1::LongAA, seq2::LongAA) = pairalign(GlobalAlignment(), seq1, seq2, score_model)

    # Calculate sequence identities from pairwise alignments
    get_denominator(a::T) where {T<:PairwiseAlignmentResult} = (a |> alignment |> length) # Alignment length
    #get_denominator(a::T) where {T<:PairwiseAlignmentResult} = minimum([length(a.aln.a.seq), length(a.aln.b)]) # Length of shorter sequence
    get_sequence_identity(a::T) where {T<:PairwiseAlignmentResult} = (a |> alignment |> count_matches) / get_denominator(a)
    #get_sequence_similarit(a::T) where {T<:PairwiseAlignmentResult} = (a)
    get_sequence_similarity(a::T) where {T<:PairwiseAlignmentResult} = count_positives(a, SUBSTITUTION_MATRIX) / get_denominator(a)

    max_sequence_identity = zeros(Float64, length(generated))
    max_sequence_similarity = zeros(Float64, length(generated))
    max_sequence_identity_index = zeros(Int, length(generated))
    max_sequence_similarity_index = zeros(Int, length(generated))
    for (i, generated_seq) in enumerate(generated)
        for (j, train_seq) in enumerate(train)
            ((self_comparison) && (i == j)) && continue # Skip self-comparison

            sequence_alignment = align_sequences(generated_seq, train_seq)
            sequence_identity = get_sequence_identity(sequence_alignment)
            if sequence_identity > max_sequence_identity[i]
                max_sequence_identity[i] = sequence_identity
                max_sequence_identity_index[i] = j
            end

            sequence_similarity = get_sequence_similarity(sequence_alignment)
            if sequence_similarity > max_sequence_similarity[i]
                max_sequence_similarity[i] = sequence_similarity
                max_sequence_similarity_index[i] = j
            end
        end
    end

    return max_sequence_identity, max_sequence_similarity, max_sequence_identity_index, max_sequence_similarity_index
end

#= function get_ncbi_blast_sequence_identity(train_path, generated_path)
    ret = read(pipeline(`blastp -query $train_path -subject $generated_path -outfmt "10 pident ppos ppos"`), String) # TODO use train, generated
    ret = map(row -> split(row, ','), split(ret, "\n")) # split rows by "\n"
    ret = map(row -> map(string_value -> replace(string_value, "\r" => ""), row), ret) # remove "\r"
    if ret[end] == [""]
        ret = ret[1:end-1]
    end
    identity = map(row -> parse(Float64, row[1]), ret)
    similarity = map(row -> parse(Float64, row[2]), ret)
    return identity, similarity
end =#

function _get_save_path(data_path::String;
	save_suffix::Union{String, Nothing} = nothing,
)
    extension = split(data_path, ".")[end]
	base_path = data_path[1:end-length(extension)-1]
	suffix = isnothing(save_suffix) ? "" : "_" * save_suffix
	save_path = base_path * "_max_sequence_identity" * suffix * ".csv"
	return save_path
end

########## Main ##########
function max_sequence_identity(
    fasta_path::String;
    train_path::Union{String,Nothing}=nothing,
)
    if isnothing(train_path)
        main_train_sequences(fasta_path)
    else
        main_generated_sequences(fasta_path=fasta_path, train_path=train_path)
    end
end

########## Main with reference sequences ##########
function main_generated_sequences(;
    fasta_path::String,
    train_path::String,
)
    generated_identifiers, generated = separate_identifiers(load_fasta_sequences(fasta_path; load_identifiers=true))
    train_identifiers, train = separate_identifiers(load_fasta_sequences(train_path; load_identifiers=true))
    save_path = _get_save_path(fasta_path)
    evaluate_max_sequence_identity(train, generated, generated_identifiers, train_identifiers; save_path)
end

########## Main without reference sequences ##########
function main_train_sequences(train_path::String)
    train_identifiers, train = separate_identifiers(load_fasta_sequences(train_path; load_identifiers=true))
    _main_train_sequences(train_path, train, train_identifiers)
end
function _main_train_sequences(train_path::String, train, train_identifiers)
    save_path = _get_save_path(train_path; save_suffix="self")
    evaluate_max_sequence_identity(train, train, train_identifiers, train_identifiers; save_path, return_second_max=true)
end
