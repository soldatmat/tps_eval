from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Sequence-branch NUMERIC targets (present for every dataset → comparison plots)
# ---------------------------------------------------------------------------
# The original hardcoded set, plus the newer sequence metrics. Targets listed
# here that lack an explicit MIN_VAL/MAX_VAL/TICKS entry are auto-ranged from
# the data at plot time (see boxplot_comparison / density_comparison).
TARGETS = [
    "sequence_identity",
    "sequence_identity_self",
    "sequence_similarity",
    "sequence_similarity_self",
    "min_embedding_distance",
    "min_embedding_distance_self",
    # --- local (MMseqs2) sequence search: gen-vs-train + within-set (_self) ---
    # identity/similarity in [0, 100]; coverage in [0, 100]. Comparative (gen vs train)
    # -> still EXCLUDED from reference bands. similarity is NaN for the MMseqs2 backend.
    "local_sequence_identity",
    "local_sequence_identity_self",
    "local_sequence_similarity",
    "local_sequence_similarity_self",
    "local_coverage",
    "local_coverage_self",
    "isTPS",
    "isTPS_seq",
    "soluble",
    # --- motif pair distance ---
    "motif_start_distance",
    "motif_gap",
    # --- ESM pseudo-perplexity (naturalness) ---
    "esm_pseudo_perplexity",
    "esm_mean_pll",
    # --- DIAMOND Swiss-Prot homology search ---
    "swissprot_top_pident",
    "swissprot_best_nontps_pident",
    "swissprot_n_tps_in_topN",
    # --- k-NN coarse-label transfer (gen-only; <gen>_knn_label_transfer.csv) ---
    # `confidence` is the ensembled calibrated confidence in [0, 1]; the predicted
    # label itself is categorical (see CATEGORICAL_TARGETS below).
    "confidence",
    # --- substrate-class combiner (gen-only; <gen>_substrate_class.csv) ---
    # `substrate_confidence` is the substrate k-NN calibrated confidence in [0, 1]
    # (renamed from the CSV's `confidence` to avoid colliding with the k-NN target);
    # `n_signals_agree` counts how many of {pocket-band, EE} corroborate the call.
    # The predicted substrate + agreement flags are categorical (see below).
    "substrate_confidence",
    "n_signals_agree",
]


LOAD = {
    "sequence_identity": ["max_sequence_identity"],
    "sequence_identity_self": ["max_sequence_identity_self"],
    "sequence_similarity": ["max_sequence_identity"],
    "sequence_similarity_self": ["max_sequence_identity_self"],
    "min_embedding_distance": ["min_embedding_distance"],
    "min_embedding_distance_self": ["min_embedding_distance_self"],
    "local_sequence_identity": ["local_sequence_search"],
    "local_sequence_similarity": ["local_sequence_search"],
    "local_coverage": ["local_sequence_search"],
    "local_sequence_identity_self": ["local_sequence_search_self"],
    "local_sequence_similarity_self": ["local_sequence_search_self"],
    "local_coverage_self": ["local_sequence_search_self"],
    "isTPS": ["enzyme_explorer"],
    "isTPS_seq": ["enzyme_explorer_sequence_only"],
    "soluble": ["soluprot"],
    "motif_start_distance": ["motif_pair_distance"],
    "motif_gap": ["motif_pair_distance"],
    "esm_pseudo_perplexity": ["esm_pseudo_perplexity"],
    "esm_mean_pll": ["esm_pseudo_perplexity"],
    "swissprot_top_pident": ["swissprot_search"],
    "swissprot_best_nontps_pident": ["swissprot_search"],
    "swissprot_n_tps_in_topN": ["swissprot_search"],
    "confidence": ["knn_label_transfer"],
    "substrate_confidence": ["substrate_class"],
    "n_signals_agree": ["substrate_class"],
}


# Axis bounds / ticks are defined only for the targets with a meaningful fixed
# scale (probabilities, identities). Anything absent here is auto-ranged.
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
    "confidence": 0.0 - 0.01,
    "substrate_confidence": 0.0 - 0.01,
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
    "confidence": 1.0 + 0.01,
    "substrate_confidence": 1.0 + 0.01,
}


# local_sequence_search identity/similarity/coverage are PERCENTAGES in [0, 100]
# (unlike the global max_sequence_identity which is a [0, 1] fraction). similarity is
# NaN for the MMseqs2 backend (no positives field) -> auto-handled by the plot layer.
for _t in (
    "local_sequence_identity", "local_sequence_identity_self",
    "local_sequence_similarity", "local_sequence_similarity_self",
    "local_coverage", "local_coverage_self",
):
    MIN_VAL[_t] = 0.0 - 1.0
    MAX_VAL[_t] = 100.0 + 1.0


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
    "confidence": _ticks(0.0, 1.0, 0.05),
    "substrate_confidence": _ticks(0.0, 1.0, 0.05),
}

# Percentage-scale ticks for the local_sequence_search targets ([0, 100]).
for _t in (
    "local_sequence_identity", "local_sequence_identity_self",
    "local_sequence_similarity", "local_sequence_similarity_self",
    "local_coverage", "local_coverage_self",
):
    TICKS[_t] = _ticks(0.0, 100.0, 10.0)


# `None` (or absent) → no threshold line drawn.
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


# Ridge-plot dataset spacing. Absent → default OFFSET_DEFAULT.
OFFSET_DEFAULT = 3.0
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


# ---------------------------------------------------------------------------
# Sequence-branch CATEGORICAL / BOOLEAN targets → grouped count plots
# ---------------------------------------------------------------------------
# `swissprot_top_is_tps` is boolean. Motif-presence columns (from the *_motifs.csv,
# whose column names ARE the regex patterns) are discovered dynamically at plot
# time, so they aren't enumerated here.
CATEGORICAL_TARGETS = {
    "swissprot_top_is_tps": ["swissprot_search"],
    # k-NN ensembled predicted coarse label (gen-only; "unknown" = abstained).
    "predicted_label": ["knn_label_transfer"],
    # Substrate-class combiner (gen-only): the fused substrate call + whether the
    # pocket-volume size band corroborates it.
    "predicted_substrate": ["substrate_class"],
    "substrate_agreement": ["substrate_class"],
}


# ---------------------------------------------------------------------------
# STRUCTURE-branch targets (structures exist for the generated set only →
# single-distribution plots, no train series). Each maps the target column to
# the `<structs_dir>_<suffix>.csv` file it is read from. Files are discovered by
# suffix in the input directory at plot time; a missing file is skipped cleanly.
# ---------------------------------------------------------------------------
# suffix -> list of NUMERIC metric columns in that CSV
STRUCTURE_NUMERIC = {
    "_plddt.csv": [
        "mean_plddt",
        "median_plddt",
        "min_plddt",
        "frac_plddt_confident",
    ],
    "_structural_identity.csv": [
        "structural_tmscore_to_known",
        "structural_lddt_to_known",
    ],
    "_motif_structural_distance.csv": [
        "motif_centroid_distance",
        "motif_min_ca_distance",
    ],
    "_active_site_geometry.csv": [
        "carboxylate_convergence_radius",
        "n_coordinating_oxygens",
        "metal_point_void",
        "catalytic_constellation_rmsd",
    ],
    "_aromatic_lining.csv": [
        "n_pocket_aromatics",
        "n_inward_facing_aromatics",
        "aromatic_fraction",
    ],
    "_diphosphate_sensor.csv": [
        "n_diphosphate_basic_residues",
        "n_RY_pairs",
    ],
    # Ion-placement check (AF3 holo folds): how well the co-folded Mg/Mn ions land in
    # the carboxylate cage. NaN / 0 for apo/ESMFold (no ions) -> a not-applicable row.
    "_ion_site_check.csv": [
        "min_ion_to_cage_dist",
        "n_ions_in_site",
        "max_coordinating_contacts",
        "n_motif_coord_asp",
        "n_motif_coord_nse",
        "mg_to_motif_dist",
    ],
    "_substrate_positioning.csv": [
        "diphosphate_to_cage_dist",
        "min_diphosphate_to_cage_oxygen",
        "diphosphate_to_nearest_ion",
        "diphosphate_to_ion_centroid",
        "reactive_carbon_to_cage_dist",
        "reactive_carbon_to_nearest_ion",
        "reactive_carbon_to_ion_centroid",
        "substrate_plddt",
    ],
    "_cyclization_geometry.csv": [
        "substrate_rgyr",
        "foldback_c1_to_distal",
        "substrate_endtoend",
        "n_aromatic_carbon_contacts",
        "frac_aromatic_track",
        "n_aromatics_lining",
        "mean_carbon_to_aromatic",
    ],
    "_radius_of_gyration.csv": [
        "radius_of_gyration",
        "asphericity",
        "acylindricity",
    ],
    "_pocket_descriptors.csv": [
        "catalytic_pocket_volume",
        "pocket_hydrophobicity",
        "pocket_enclosure",
        "pocket_n_alpha_spheres",
        "pocket_depth",
        "pocket_sasa_per_volume",
        "p2rank_catalytic_site_score",
        "p2rank_catalytic_pocket_rank",
    ],
    "_domain_composition.csv": [
        "n_domains",
        "alpha",
        "beta",
        "gamma",
        "delta",
        "epsilon",
        "ids",
        "terpene_synth_C",
    ],
    "_aggregation.csv": [
        "a3d_avg_score",
        "a3d_total_pos_score",
        "a3d_max_score",
    ],
    "_foldseek_swissprot_search.csv": [
        "foldseek_sprot_top_tmscore",
        "foldseek_sprot_best_nontps_tmscore",
        "foldseek_sprot_n_tps_in_topN",
    ],
    "_proteinmpnn_score.csv": [
        "proteinmpnn_nll",
    ],
    "_global_confidence.csv": [
        "ptm",
        "iptm",
    ],
    "_interdomain_pae.csv": [
        "mean_interdomain_pae",
        "max_interdomain_pae",
    ],
    # Opt-in; the CSV is frequently absent (skipped cleanly).
    "_self_consistency.csv": [
        "sc_rmsd_min",
        "sc_rmsd_mean",
    ],
    # SDR specificity-divergence: global similarity to the nearest known TPS and the
    # fraction-identical fraction at the specificity-determining residues. Both in
    # [0, 1] (fixed-scale below); the divergence flag itself is categorical.
    "_sdr_divergence.csv": [
        "nearest_neighbour_similarity",
        "sdr_identity",
    ],
    # Domain-level structural identity: best TM-score / lddt of each design's detected
    # TPS domains to the known martsDB reference domains, plus how many domains EE found.
    "_domain_structural_identity.csv": [
        "domain_structural_tmscore_to_known",
        "domain_structural_lddt_to_known",
        "n_detected_domains",
    ],
}

# suffix -> list of CATEGORICAL / BOOLEAN metric columns in that CSV
STRUCTURE_CATEGORICAL = {
    "_domain_composition.csv": ["domain_architecture"],
    "_foldseek_swissprot_search.csv": ["foldseek_sprot_top_is_tps"],
    "_diphosphate_sensor.csv": ["has_RY_pair"],
    "_ion_site_check.csv": ["ion_in_site", "well_placed", "mg_canonical_motif_coordination"],
    "_substrate_positioning.csv": ["substrate_present", "substrate_in_site"],
    "_cyclization_geometry.csv": ["substrate_present"],
    "_sdr_divergence.csv": ["specificity_divergence"],
}


# Fixed-scale axis bounds for the few structure metrics that have a natural
# [0, 1] domain (everything else is auto-ranged).
_PROB_LIKE_STRUCTURE = (
    "frac_plddt_confident",
    "structural_tmscore_to_known",
    "structural_lddt_to_known",
    "foldseek_sprot_top_tmscore",
    "foldseek_sprot_best_nontps_tmscore",
    "nearest_neighbour_similarity",
    "sdr_identity",
    "domain_structural_tmscore_to_known",
    "domain_structural_lddt_to_known",
)
for _t in _PROB_LIKE_STRUCTURE:
    MIN_VAL[_t] = 0.0 - 0.01
    MAX_VAL[_t] = 1.0 + 0.01
    TICKS[_t] = _ticks(0.0, 1.0, 0.05)
