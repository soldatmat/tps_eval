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

LOAD = Dict([
    ("sequence_identity", [:max_sequence_identity]),
    ("sequence_identity_self", [:max_sequence_identity_self]),
    ("sequence_similarity", [:max_sequence_identity]),
    ("sequence_similarity_self", [:max_sequence_identity_self]),
    ("min_embedding_distance", [:min_embedding_distance]),
    ("min_embedding_distance_self", [:min_embedding_distance_self]),
    ("isTPS", [:enzyme_explorer]),
    ("isTPS_seq", [:enzyme_explorer_sequence_only]),
    ("soluble", [:soluprot]),
])

MIN_VAL = Dict([
    ("sequence_identity", 0.0 - 0.01),
    ("sequence_identity_self", 0.0 - 0.01),
    ("sequence_similarity", 0.0 - 0.01),
    ("sequence_similarity_self", 0.0 - 0.01),
    ("min_embedding_distance", 0.0 - 0.05),
    ("min_embedding_distance_self", 0.0 - 0.05),
    ("isTPS", 0.0 - 0.01),
    ("isTPS_seq", 0.0 - 0.01),
    ("soluble", 0.0 - 0.01),
])


MAX_VAL = Dict([
    ("sequence_identity", 1.0 + 0.01),
    ("sequence_identity_self", 1.0 + 0.01),
    ("sequence_similarity", 1.0 + 0.01),
    ("sequence_similarity_self", 1.0 + 0.01),
    ("min_embedding_distance", 6.0),
    ("min_embedding_distance_self", 6.0),
    ("isTPS", 1.0 + 0.01),
    ("isTPS_seq", 1.0 + 0.01),
    ("soluble", 1.0 + 0.01),
])

TICKS = Dict([
    ("sequence_identity", 0.0:0.05:1.0),
    ("sequence_identity_self", 0.0:0.05:1.0),
    ("sequence_similarity", 0.0:0.05:1.0),
    ("sequence_similarity_self", 0.0:0.05:1.0),
    ("min_embedding_distance", 0.0:0.25:6.0),
    ("min_embedding_distance_self", 0.0:0.25:6.0),
    ("isTPS", 0.0:0.05:1.0),
    ("isTPS_seq", 0.0:0.05:1.0),
    ("soluble", 0.0:0.05:1.0),
])

THRESHOLD = Dict([
    ("sequence_identity", 0.5),
    ("sequence_identity_self", nothing),
    ("sequence_similarity", 0.55),
    ("sequence_similarity_self", nothing),
    ("min_embedding_distance", 1.25),
    ("min_embedding_distance_self", nothing),
    ("isTPS", 0.35),
    ("isTPS_seq", 0.35),
    ("soluble", 0.5),
])

OFFSET = Dict([ # TODO fill values
    ("sequence_identity", 3.0),
    ("sequence_identity_self", 3.0),
    ("sequence_similarity", 3.0),
    ("sequence_similarity_self", 3.0),
    ("min_embedding_distance", 3.0),
    ("min_embedding_distance_self", 3.0),
    ("isTPS", 3.0),
    ("isTPS_seq", 3.0),
    ("soluble", 3.0),
])
