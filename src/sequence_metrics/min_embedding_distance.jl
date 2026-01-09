using Distances

include("../data/embeddings.jl")

# Vector{Vector{Union{Float64, Missing}}} -> Vector{Vector{Float64}}
preprocess_embeddings(embeddings_df) = Vector{Vector{typeof(embeddings_df.embedding[1][1])}}(embeddings_df.embedding)
get_distances(embeddings1, embeddings2) = pairwise(Euclidean(), embeddings1, embeddings2)
get_min_distances(distances) = map(dists -> findmin(dists), eachrow(distances))
get_max_distances(distances) = map(dists -> findmax(dists), eachrow(distances))

function save_embeddings(
	ids::Vector{<:AbstractString},
	min_dist::Vector{Float64},
	min_dist_hits::Vector{<:AbstractString},
	save_path::String,
)
	df = DataFrame(ID = ids, min_embedding_distance = min_dist, min_embedding_distance_hit = min_dist_hits)
	CSV.write(save_path, df)
end

function _min_embedding_distance(
	train_df,
	generated_df;
	save_path::Union{String, Nothing} = nothing,
)
	generated_embeddings = preprocess_embeddings(generated_df)
    train_embeddings = preprocess_embeddings(train_df)

	distances = get_distances(generated_embeddings, train_embeddings)

	min_dist_results = get_min_distances(distances)
	min_dist = [res[1] for res in min_dist_results]
	min_dist_index = [res[2] for res in min_dist_results]
	min_dist_hits = train_df.ID[min_dist_index]

	isnothing(save_path) || save_embeddings(
		generated_df.ID,
		min_dist,
		min_dist_hits,
		save_path,
	)

	return min_dist
end


function _min_embedding_distance_self(
	train_df;
	save_path::Union{String, Nothing} = nothing,
)
    train_embeddings = preprocess_embeddings(train_df)

	distances = get_distances(train_embeddings, train_embeddings)
	# Set diagonal to Inf to ignore self-distance
	for i in 1:size(distances, 1)
		distances[i, i] = Inf
	end

	min_dist, min_dist_index = get_min_distances(distances)
	min_dist_hits = train_df.ID[min_dist_index]

	isnothing(save_path) || save_embeddings(
		train_df.ID,
		min_dist,
		min_dist_hits,
		save_path,
	)

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
	_min_embedding_distance_self(train_df;
		save_path,
	)
end
