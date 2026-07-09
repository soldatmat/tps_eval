using FASTX

#=
for record in reader
    record |> identifier
    record |> sequence
end
=#

function load_fasta_sequences(file_path::String; remove_padding::Bool=true, load_identifiers::Bool=false)
    reader = FASTAReader(open(file_path))
    sequences = read_sequences(reader; remove_padding, load_identifiers)
    close(reader)
    return sequences
end

function separate_identifiers(train)
    train_identifiers = getfield.(train, 1)
    train = getfield.(train, 2)
    return train_identifiers, train
end

function load_fasta_msa(file_path::String; load_identifiers::Bool=false)
    reader = FASTAReader(open(file_path))
    msa = read_sequences(reader; remove_padding=false, load_identifiers)
    close(reader)
    reduce(hcat, collect.(msa)) |> permutedims
end

function read_sequences(reader::FASTX.FASTA.Reader; remove_padding::Bool=true, load_identifiers::Bool=false)
    sequences = load_identifiers ? Vector{Tuple{String,String}}([]) : Vector{String}([])
    for record in reader
        seq = record |> FASTX.sequence |> string
        if remove_padding
            seq = replace(seq, "-" => "")
        end
        if load_identifiers
            id = record |> identifier |> string
        end
        load_identifiers ? append!(sequences, [(id, seq)]) : append!(sequences, [seq])
    end
    return sequences
end
