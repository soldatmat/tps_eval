using CSV
using DataFrames

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
