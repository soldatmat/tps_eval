using CSV
using DataFrames

include("sequences.jl")


DEFAULT_LOAD = [
    :sequence,
    :max_sequence_identity,
    :max_sequence_identity_self,
    :embedding,
    :min_embedding_distance,
    :min_embedding_distance_self,
    :enzyme_explorer_sequence_only,
    :enzyme_explorer,
    :motifs,
    :soluprot,
]

type_to_substrate = Dict([
    ("mono", ["CC(C)=CCCC(C)=CCOP([O_])(=O)OP([O_])([O_])=O", "CC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O (Geranyl pyrophosphate)"]), # GPP
    ("sesq", ["CC(C)=CCCC(C)=CCCC(C)=CCOP([O_])(=O)OP([O_])([O_])=O", "CC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O (Farnesyl pyrophosphate)"]), # FPP
    ("di", ["CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O_])(=O)OP([O_])([O_])=O", "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O (Geranylgeranyl pyrophosphate)"]), # GGPP
    ("sester", ["CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O_])(=O)OP([O_])([O_])=O", "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O (Geranylfarnesyl pyrophosphate)"]), # GFPP
    ("tri", ["CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC1OC1(C)C", "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC1OC1(C)C ((S)-2,3-epoxysqualene)"]), # Epoxysqualene (also could use Squalene)
])


function get_target_substrate(tps_type, df_names)
    target_substrate_names = type_to_substrate[tps_type]
    for name in target_substrate_names
        if name in df_names
            return name
        end
    end
    return nothing
end

function construct_result_paths(fasta_path::String)
    if endswith(fasta_path, ".fasta")
        remove_length = 6
    elseif endswith(fasta_path, ".fa")
        remove_length = 3
    else
        error("The fasta file path must end with '.fasta' or '.fa'.")
    end
    partial_sequences_path = fasta_path[1:end-remove_length]

    max_sequence_identity_path = partial_sequences_path * "_max_sequence_identity.csv"
    max_sequence_identity_self_path = partial_sequences_path * "_max_sequence_identity_self.csv"

    embedding_esm1b_path = partial_sequences_path * "_embedding_esm1b.csv"
    min_embedding_distance_esm1b_path = partial_sequences_path * "_embedding_esm1b_min_embedding_distance.csv"
    min_embedding_distance_esm1b_self_path = partial_sequences_path * "_embedding_esm1b_min_embedding_distance_self.csv"
    
    enzyme_explorer_sequence_only_path = partial_sequences_path * "_enzyme_explorer_sequence_only.csv"
    enzyme_explorer_path = partial_sequences_path * "_enzyme_explorer.csv"

    motifs_path = partial_sequences_path * "_motifs.csv"

    soluprot_path = partial_sequences_path * "_soluprot.csv"

    return (
        fasta_path,
        max_sequence_identity_path, 
        max_sequence_identity_self_path,
        embedding_esm1b_path, 
        min_embedding_distance_esm1b_path, 
        min_embedding_distance_esm1b_self_path,
        enzyme_explorer_sequence_only_path,
        enzyme_explorer_path,
        motifs_path,
        soluprot_path
    )
end

function load_results(
    fasta_path::String,
    max_sequence_identity_path::String, 
    max_sequence_identity_self_path::String,
    embedding_esm1b_path::String, 
    min_embedding_distance_esm1b_path::String, 
    min_embedding_distance_esm1b_self_path::String,
    enzyme_explorer_sequence_only_path::String,
    enzyme_explorer_path::String,
    motifs_path::String, # TODO add
    soluprot_path::String;

    tps_type::Union{String, Nothing}=nothing,
    load::Union{Vector{Symbol},Nothing}=DEFAULT_LOAD,
)
    if isnothing(load)
        load = DEFAULT_LOAD
    end

    dataframes = Vector{DataFrame}([])

    # Sequences
    if :sequence in load
        records = load_fasta_sequences(fasta_path; remove_padding=true, load_identifiers=true)
        sequences = DataFrame(map(idx -> getindex.(records, idx), eachindex(first(records))), [:ID, :sequence])
        sequences.msa = load_fasta_sequences(fasta_path; remove_padding=false)
        push!(dataframes, sequences)
    end

    # ESM1b embeddings
    if :embedding in load
        embedding = load_embeddings(embedding_esm1b_path)
        push!(dataframes, embedding)
    end

    # EnzymeExplorer predictions
    if :enzyme_explorer_sequence_only in load
        scores = CSV.read(enzyme_explorer_sequence_only_path, DataFrame)
        enzyme_explorer_column_names = ["ID", "isTPS"]
        if !isnothing(tps_type)
            target_substrate = get_target_substrate(tps_type, names(scores))
            push!(enzyme_explorer_column_names, target_substrate)
        end
        rename!(scores, Dict(col => strip(col) for col in names(scores)))
        scores = scores[!, enzyme_explorer_column_names]
        isnothing(tps_type) || rename!(scores, target_substrate => :target_substrate)
        rename!(scores, :isTPS => :isTPS_seq)
        push!(dataframes, scores)
    end

    if :enzyme_explorer in load
        scores = CSV.read(enzyme_explorer_path, DataFrame)
        enzyme_explorer_column_names = ["ID", "isTPS"]
        if !isnothing(tps_type)
            target_substrate = get_target_substrate(tps_type, names(scores))
            push!(enzyme_explorer_column_names, target_substrate)
        end
        rename!(scores, Dict(col => strip(col) for col in names(scores)))
        scores = scores[!, enzyme_explorer_column_names]
        isnothing(tps_type) || rename!(scores, target_substrate => :target_substrate)
        push!(dataframes, scores)
    end

    # Max sequence identity & similarity to train sequences
    if :max_sequence_identity in load
        max_sequence_identity = CSV.read(max_sequence_identity_path, DataFrame)
        push!(dataframes, max_sequence_identity)
    end

    # Max sequence identity & similarity to other sequences
    if :max_sequence_identity_self in load
        max_sequence_identity_self = CSV.read(max_sequence_identity_self_path, DataFrame)
        rename!(max_sequence_identity_self, :sequence_identity => :sequence_identity_self)
        rename!(max_sequence_identity_self, :sequence_similarity => :sequence_similarity_self)
        push!(dataframes, max_sequence_identity_self)
    end

    # Min ESM1b embedding distances to train sequences
    if :min_embedding_distance in load
        min_embedding_distance = CSV.read(min_embedding_distance_esm1b_path, DataFrame)
        push!(dataframes, min_embedding_distance)
    end

    # Min ESM1b embedding distances to other sequences
    if :min_embedding_distance_self in load
        min_embedding_distance_self = CSV.read(min_embedding_distance_esm1b_self_path, DataFrame)
        rename!(min_embedding_distance_self, :min_embedding_distance => :min_embedding_distance_self)
        push!(dataframes, min_embedding_distance_self)
    end

    # Presence of motifs
    if :motifs in load
        motifs = CSV.read(motifs_path, DataFrame)
        select!(motifs, Not(:sequence))
        push!(dataframes, motifs)
    end

    # SoluProt solubility predictions
    if :soluprot in load
        solubility = CSV.read(soluprot_path, DataFrame)
        rename!(solubility, :fa_id => :ID)
        push!(dataframes, solubility)
    end

    # Join dataframes
    # TODO problem with whitespaces in ID (e.g. "SRA.SRR18681978 (AgTS-1)")
    # some of the files include the full id ("SRA.SRR18681978 (AgTS-1)") and some strip it ("SRA.SRR18681978")
    
    # Keep only the first part of the ID (e.g.
    #    "SRA.SRR18681978" instead of "SRA.SRR18681978 (AgTS-1)",
    #    "seq7235 instead of "seq7235 traget 10"
    # )
    for df in dataframes
        if "ID" in names(df)
            println("Warning: stripping ID column of parts after whitespaces.")
            df.ID .= split.(df.ID, ' ').|>first
        end
    end

    df = length(dataframes) == 1 ? dataframes[1] : outerjoin(dataframes...; on=:ID)

    return df
end

function load_embeddings(file_path::String)
    original_df = CSV.read(file_path, DataFrame)
    sequence_embeddings = [collect(values(row)) for row in eachrow(original_df[!, Not(:id)])]

    embedding_df = DataFrame(
        id = original_df.id,
        embedding = sequence_embeddings,
    )
    rename!(embedding_df, :id => :ID)

    return embedding_df
end

function load_results(
    fasta_path::String;
    tps_type::Union{String, Nothing}=nothing,
    load::Union{Vector{Symbol},Nothing}=DEFAULT_LOAD,
)
    (
        fasta_path,
        max_sequence_identity_path, 
        max_sequence_identity_self_path,
        embedding_esm1b_path, 
        min_embedding_distance_esm1b_path, 
        min_embedding_distance_esm1b_self_path,
        enzyme_explorer_sequence_only_path,
        enzyme_explorer_path,
        motifs_path,
        soluprot_path
    ) = construct_result_paths(fasta_path)

    return load_results(
        fasta_path,
        max_sequence_identity_path, 
        max_sequence_identity_self_path,
        embedding_esm1b_path, 
        min_embedding_distance_esm1b_path, 
        min_embedding_distance_esm1b_self_path,
        enzyme_explorer_sequence_only_path,
        enzyme_explorer_path,
        motifs_path,
        soluprot_path;
        tps_type=tps_type,
        load=load,
    )
end
