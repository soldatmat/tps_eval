from __future__ import annotations

import numpy as np


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


LOAD = {
    "sequence_identity": ["max_sequence_identity"],
    "sequence_identity_self": ["max_sequence_identity_self"],
    "sequence_similarity": ["max_sequence_identity"],
    "sequence_similarity_self": ["max_sequence_identity_self"],
    "min_embedding_distance": ["min_embedding_distance"],
    "min_embedding_distance_self": ["min_embedding_distance_self"],
    "isTPS": ["enzyme_explorer"],
    "isTPS_seq": ["enzyme_explorer_sequence_only"],
    "soluble": ["soluprot"],
}


MIN_VAL = {
    "sequence_identity": 0.0 - 0.01,
    "sequence_identity_self": 0.0 - 0.01,
    "sequence_similarity": 0.0 - 0.01,
    "sequence_similarity_self": 0.0 - 0.01,
    "min_embedding_distance": 0.0 - 0.05,
    "min_embedding_distance_self": 0.0 - 0.05,
    "isTPS": 0.0 - 0.01,
    "isTPS_seq": 0.0 - 0.01,
    "soluble": 0.0 - 0.01,
}


MAX_VAL = {
    "sequence_identity": 1.0 + 0.01,
    "sequence_identity_self": 1.0 + 0.01,
    "sequence_similarity": 1.0 + 0.01,
    "sequence_similarity_self": 1.0 + 0.01,
    "min_embedding_distance": 6.0,
    "min_embedding_distance_self": 6.0,
    "isTPS": 1.0 + 0.01,
    "isTPS_seq": 1.0 + 0.01,
    "soluble": 1.0 + 0.01,
}


def _ticks(start: float, stop: float, step: float) -> np.ndarray:
    return np.round(np.arange(start, stop + step / 2, step), 10)


TICKS = {
    "sequence_identity": _ticks(0.0, 1.0, 0.05),
    "sequence_identity_self": _ticks(0.0, 1.0, 0.05),
    "sequence_similarity": _ticks(0.0, 1.0, 0.05),
    "sequence_similarity_self": _ticks(0.0, 1.0, 0.05),
    "min_embedding_distance": _ticks(0.0, 6.0, 0.25),
    "min_embedding_distance_self": _ticks(0.0, 6.0, 0.25),
    "isTPS": _ticks(0.0, 1.0, 0.05),
    "isTPS_seq": _ticks(0.0, 1.0, 0.05),
    "soluble": _ticks(0.0, 1.0, 0.05),
}


THRESHOLD = {
    "sequence_identity": 0.5,
    "sequence_identity_self": None,
    "sequence_similarity": 0.55,
    "sequence_similarity_self": None,
    "min_embedding_distance": 1.25,
    "min_embedding_distance_self": None,
    "isTPS": 0.35,
    "isTPS_seq": 0.35,
    "soluble": 0.5,
}


OFFSET = {
    "sequence_identity": 3.0,
    "sequence_identity_self": 3.0,
    "sequence_similarity": 3.0,
    "sequence_similarity_self": 3.0,
    "min_embedding_distance": 3.0,
    "min_embedding_distance_self": 3.0,
    "isTPS": 3.0,
    "isTPS_seq": 3.0,
    "soluble": 3.0,
}
