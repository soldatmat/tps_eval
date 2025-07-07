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
    generated_identifiers;
    save_path::Union{String,Nothing}=nothing,
    return_second_max::Bool=false
)
    # --- Maximum sequence identity to any train sequence --------------------------------------------------
    sequence_identity, sequence_similarity = get_pairwise_sequence_identity_with_pairwise_alignment(train, generated)

    # Get maximum sequence identity per generated sequence
    get_max_identity(sequence_identity) = map(identities -> maximum(identities), eachcol(sequence_identity))
    get_second_max_identity(sequence_identity) = getindex.(sort.(eachcol(sequence_identity), rev=true), 2)

    max_sequence_identity = return_second_max ? get_second_max_identity(sequence_identity) : get_max_identity(sequence_identity)
    max_sequence_similarity = return_second_max ? get_second_max_identity(sequence_similarity) : get_max_identity(sequence_similarity)

    # Save calculated novelty into CSV
    if !isnothing(save_path)
        df = DataFrame(
            :ID => generated_identifiers,
            :sequence_identity => max_sequence_identity,
            :sequence_similarity => max_sequence_similarity,
        )
        CSV.write(save_path, df)
    end

    return generated_identifiers, max_sequence_identity, max_sequence_similarity
end

"""
Computes sequence identity between each pair of sequences from `train` and `generated`

Sequence identity implementation:
- constructs pairwise alignment for each pair of sequences
- identity = #Matches / length(alignment)
"""
function get_pairwise_sequence_identity_with_pairwise_alignment(train::AbstractVector{LongAA}, generated::AbstractVector{LongAA})
    # Construct pairwise alignments of each generated sequence to each train sequence
    score_model = AffineGapScoreModel(SUBSTITUTION_MATRIX, gap_open=-11, gap_extend=-1)
    align_sequences(seq1::LongAA, seq2::LongAA) = pairalign(GlobalAlignment(), seq1, seq2, score_model)
    sequence_alignments = mapreduce(generated_seq -> align_sequences.(Ref(generated_seq), train), hcat, generated)
    # `sequence_alignments` is length(train)×length(generated) Matrix{PairwiseAlignmentResult{Int64, LongAA, LongAA}}
    # -> sequence_alignments[:,1] are pairwise alignments of 1st generated sequence to each train sequence

    # Calculate sequence identities from pairwise alignments
    get_denominator(a::T) where {T<:PairwiseAlignmentResult} = (a |> alignment |> length) # Alignment length
    #get_denominator(a::T) where {T<:PairwiseAlignmentResult} = minimum([length(a.aln.a.seq), length(a.aln.b)]) # Length of shorter sequence
    get_sequence_identity(a::T) where {T<:PairwiseAlignmentResult} = (a |> alignment |> count_matches) / get_denominator(a)
    #get_sequence_similarit(a::T) where {T<:PairwiseAlignmentResult} = (a)
    get_sequence_similarity(a::T) where {T<:PairwiseAlignmentResult} = count_positives(a, SUBSTITUTION_MATRIX) / get_denominator(a)

    sequence_identity = get_sequence_identity.(sequence_alignments)
    sequence_similarity = get_sequence_similarity.(sequence_alignments)

    return sequence_identity, sequence_similarity
end

function get_pairwise_sequence_identity_with_pairwise_alignment(train, generated)
    # Convert sequences to ::LongAA
    train = LongAA.(train)
    generated = LongAA.(generated)
    get_pairwise_sequence_identity_with_pairwise_alignment(train, generated)
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
    evaluate_max_sequence_identity(train, generated, generated_identifiers; save_path)
end

########## Main without reference sequences ##########
function main_train_sequences(train_path::String)
    train_identifiers, train = separate_identifiers(load_fasta_sequences(train_path; load_identifiers=true))
    _main_train_sequences(train_path, train, train_identifiers)
end
function _main_train_sequences(train_path::String, train, train_identifiers)
    save_path = _get_save_path(train_path; save_suffix="self")
    evaluate_max_sequence_identity(train, train, train_identifiers; save_path, return_second_max=true)
end
