using Distances

include("../data/embeddings.jl")

get_min_distances(distances) = map(dists -> minimum(dists), eachrow(distances))
get_second_min_distances(distances) = getindex.(sort.(eachrow(distances), rev = false), 2)
get_max_distances(distances) = map(dists -> maximum(dists), eachrow(distances))
get_second_max_distances(distances) = getindex.(sort.(eachrow(distances), rev = true), 2)

function _min_embedding_distance(
	train_df,
	generated_df;
	save_path::Union{String, Nothing} = nothing,
	return_second_min::Bool = false,
)
	# Vector{Vector{Union{Float64, Missing}}} -> Vector{Vector{Float64}}
	generated_embeddings = Vector{Vector{typeof(generated_df.embedding[1][1])}}(generated_df.embedding)
    train_embeddings = Vector{Vector{typeof(train_df.embedding[1][1])}}(train_df.embedding)

	distances = pairwise(Euclidean(), generated_embeddings, train_embeddings)

	min_dist = return_second_min ? get_second_min_distances(distances) : get_min_distances(distances)

	# Save calculated novelty into CSV
	df = DataFrame(ID = generated_df.ID, min_embedding_distance = min_dist)
	isnothing(save_path) || CSV.write(save_path, df)

	return min_dist
end

function _get_save_path(embeddings_path::String;
	save_suffix::Union{String, Nothing} = nothing,
)
	extension = split(embeddings_path, ".")[end]
	base_path = embeddings_path[1:end-length(extension)-1]
	suffix = isnothing(save_suffix) ? "" : "_" * save_suffix
	save_path = base_path * "_min_embedding_distance" * suffix * ".csv"
	return save_path
end

########## Main ##########
function min_embedding_distance(
	embeddings_path::String;
    train_embeddings_path::Union{String,Nothing}=nothing,
	save::Bool = true,
)
	if isnothing(train_embeddings_path)
		main_train_sequences(embeddings_path; save = save)
	else
		main_generated_sequences(embeddings_path, train_embeddings_path; save = save)
	end
end

########## Main with reference sequences ##########
function main_generated_sequences(generated_embeddings_path::String, train_embeddings_path::String;
	save::Bool = true,
)
	train_df = load_embeddings(train_embeddings_path)
	generated_df = load_embeddings(generated_embeddings_path)

	save_path = _get_save_path(generated_embeddings_path)
	_min_embedding_distance(train_df, generated_df;
		save_path = save ? save_path : nothing,
		return_second_min = false,
	)
end

########## Main without reference sequences ##########
function main_train_sequences(embeddings_path::String;
	save::Bool = true,
)
	train_df = load_embeddings(embeddings_path)

	save_path = _get_save_path(embeddings_path; save_suffix = "self")
	main_train_sequences(train_df;
		save_path = save ? save_path : nothing,
	)
end

function main_train_sequences(train_df;
	save_path::Union{String, Nothing} = nothing,
)
	_min_embedding_distance(train_df, train_df;
		save_path,
		return_second_min = true,
	)
end
