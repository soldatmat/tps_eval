"""Per-metric metadata for the dashboard: a one-line explanation (hover "?"),
the mathematical range of each numeric column (shown in the column header), and
the category each metric belongs to.

Compiled from docs/TOOLS.md / README.md / the band JSONs. Adding a new metric
here is optional — the dashboard still renders metrics with no entry (no tooltip,
no range chip, and they fall into the "Other" category).
"""

from __future__ import annotations

# Display order of the metric categories (left-to-right grouping in the dashboard).
CATEGORY_ORDER = ["Fold & confidence", "Active site", "Sequence", "Function", "Novelty", "Other"]

# metric -> category. Metrics absent here fall into "Other".
METRIC_CATEGORY = {
    # Fold & confidence
    "plddt": "Fold & confidence",
    "global_confidence": "Fold & confidence",
    "interdomain_pae": "Fold & confidence",
    "radius_of_gyration": "Fold & confidence",
    "domain_composition": "Fold & confidence",
    "aggregation": "Fold & confidence",
    "self_consistency": "Fold & confidence",
    # Active site
    "active_site_geometry": "Active site",
    "pocket_descriptors": "Active site",
    "aromatic_lining": "Active site",
    "diphosphate_sensor": "Active site",
    "ion_site_check": "Active site",
    "substrate_positioning": "Active site",
    "motif_pair_distance": "Active site",
    "motif_structural_distance": "Active site",
    "motif_search": "Active site",
    # Sequence
    "esm_pseudo_perplexity": "Sequence",
    "soluprot": "Sequence",
    "proteinmpnn_score": "Sequence",
    # Function
    "enzyme_explorer_sequence_only": "Function",
    "substrate_class": "Function",
    "knn_label_transfer": "Function",
    # Novelty (comparative — similarity to a reference set)
    "max_sequence_identity": "Novelty",
    "min_embedding_distance": "Novelty",
    "local_sequence_search": "Novelty",
    "swissprot_search": "Novelty",
    "foldseek_swissprot_search": "Novelty",
    "structural_identity": "Novelty",
    "domain_structural_identity": "Novelty",
    "sdr_divergence": "Novelty",
}

METRIC_INFO = {
    # ----------------------- BANDED METRICS -----------------------
    "plddt": {
        "explanation": "Per-residue AlphaFold/ESMFold fold-confidence (pLDDT); higher means the backbone is more confidently predicted.",
        "columns": {
            "mean_plddt": "0–100",
            "median_plddt": "0–100",
            "min_plddt": "0–100",
            "frac_plddt_confident": "0–1",
            "n_residues": "integer ≥ 1",
        },
    },
    "global_confidence": {
        "explanation": "Whole-fold confidence scores (pTM and, for complexes, ipTM) summarizing how reliable the overall predicted structure is.",
        "columns": {
            "ptm": "0–1",
            "iptm": "0–1",
        },
    },
    "interdomain_pae": {
        "explanation": "Predicted alignment error between domains; low values mean the relative orientation of the design's domains is confidently placed.",
        "columns": {
            "mean_interdomain_pae": "≥ 0 (Å)",
            "max_interdomain_pae": "≥ 0 (Å)",
            "n_domains": "integer ≥ 0",
        },
    },
    "esm_pseudo_perplexity": {
        "explanation": "Sequence naturalness under ESM's masked language model; lower perplexity means the sequence is more in-distribution for real proteins.",
        "columns": {
            "esm_pseudo_perplexity": "≥ 1",
            "esm_mean_pll": "≤ 0",
            "n_residues": "integer ≥ 1",
        },
    },
    "proteinmpnn_score": {
        "explanation": "ProteinMPNN negative log-likelihood of the design's own sequence given its backbone; lower means the sequence fits the fold better.",
        "columns": {
            "proteinmpnn_nll": "≥ 0",
            "proteinmpnn_score_designed": "≥ 0",
        },
    },
    "radius_of_gyration": {
        "explanation": "Global size and shape of the fold (compactness and elongation), to spot designs that are too extended or oddly shaped.",
        "columns": {
            "radius_of_gyration": "≥ 0 (Å)",
            "asphericity": "≥ 0 (Å²)",
            "acylindricity": "≥ 0 (Å²)",
            "principal_radius_1": "≥ 0 (Å)",
            "principal_radius_2": "≥ 0 (Å)",
            "principal_radius_3": "≥ 0 (Å)",
            "n_residues": "integer ≥ 1",
        },
    },
    "active_site_geometry": {
        "explanation": "Side-chain geometry of the catalytic carboxylate cage; checks whether the metal-coordinating oxygens converge on a single metal-binding locus.",
        "columns": {
            "carboxylate_convergence_radius": "≥ 0 (Å)",
            "n_coordinating_oxygens": "integer ≥ 0",
            "metal_point_void": "≥ 0 (Å)",
            "catalytic_constellation_rmsd": "≥ 0 (Å)",
            "n_residues": "integer ≥ 1",
        },
    },
    "pocket_descriptors": {
        "explanation": "Geometry of the catalytic cavity (volume, depth, enclosure); catalytic-pocket volume is the molecular-ruler signal for product chain length.",
        "columns": {
            "metal_point_found": "",
            "catalytic_pocket_volume": "≥ 0 (Å³)",
            "pocket_hydrophobicity": "unbounded",
            "pocket_enclosure": "≥ 0",
            "pocket_n_alpha_spheres": "integer ≥ 0",
            "pocket_total_sasa": "≥ 0 (Å²)",
            "pocket_depth": "≥ 0 (Å)",
            "pocket_sasa_per_volume": "≥ 0 (Å⁻¹)",
            "fpocket_catalytic_pocket_found": "",
            "p2rank_catalytic_site_score": "≥ 0",
            "p2rank_catalytic_pocket_rank": "integer ≥ 1",
            "p2rank_catalytic_pocket_found": "",
            "n_residues": "integer ≥ 1",
        },
    },
    "aromatic_lining": {
        "explanation": "Aromatic residues (Trp/Tyr/Phe) lining the catalytic pocket that stabilize carbocation intermediates; a proxy for cyclization capability.",
        "columns": {
            "metal_point_found": "",
            "n_pocket_residues": "integer ≥ 0",
            "n_pocket_aromatics": "integer ≥ 0",
            "n_trp": "integer ≥ 0",
            "n_tyr": "integer ≥ 0",
            "n_phe": "integer ≥ 0",
            "n_his": "integer ≥ 0",
            "aromatic_fraction": "0–1",
            "n_inward_facing_aromatics": "integer ≥ 0",
            "n_residues": "integer ≥ 1",
        },
    },
    "diphosphate_sensor": {
        "explanation": "Basic residues (Arg/Lys) and the conserved RY pair that anchor and ionize the substrate diphosphate at the metal site.",
        "columns": {
            "metal_point_found": "",
            "n_diphosphate_basic_residues": "integer ≥ 0",
            "n_arg": "integer ≥ 0",
            "n_lys": "integer ≥ 0",
            "has_RY_pair": "",
            "n_RY_pairs": "integer ≥ 0",
            "n_residues": "integer ≥ 1",
        },
    },
    "motif_pair_distance": {
        "explanation": "Residue separation along the sequence between the two metal-binding motifs; a coarse check that active-site spacing is plausible.",
        "columns": {
            "motif_start_distance": "integer",
            "motif_gap": "integer",
        },
    },
    "motif_structural_distance": {
        "explanation": "3D distance between the two metal-binding motifs, approximating the span of the active-site metal cluster.",
        "columns": {
            "motif_centroid_distance": "≥ 0 (Å)",
            "motif_min_ca_distance": "≥ 0 (Å)",
            "n_residues": "integer ≥ 1",
        },
    },
    "domain_composition": {
        "explanation": "Count and type of TPS structural domains detected in each design (alpha/beta/gamma/delta/epsilon/zeta/IDS).",
        "columns": {
            "n_domains": "integer ≥ 0",
            "n_alpha": "integer ≥ 0",
            "n_beta": "integer ≥ 0",
            "n_gamma": "integer ≥ 0",
            "n_ids": "integer ≥ 0",
            "n_delta": "integer ≥ 0",
            "n_epsilon": "integer ≥ 0",
            "n_zeta": "integer ≥ 0",
            "domain_architecture": "",
        },
    },
    "aggregation": {
        "explanation": "Aggrescan3D structure-based aggregation propensity; positive scores flag surface-exposed hydrophobic patches likely to aggregate.",
        "columns": {
            "a3d_avg_score": "unbounded",
            "a3d_total_score": "unbounded",
            "a3d_max_score": "unbounded",
            "a3d_min_score": "unbounded",
            "a3d_total_pos_score": "≥ 0",
            "n_residues": "integer ≥ 1",
        },
    },
    "ion_site_check": {
        "explanation": "Whether AlphaFold3 co-folded Mg/Mn ions actually land in and are coordinated by the catalytic carboxylate cage (AF3 holo folds only).",
        "columns": {
            "metal_point_found": "",
            "n_ions_modelled": "integer ≥ 0",
            "min_ion_to_cage_dist": "≥ 0 (Å)",
            "n_ions_in_site": "integer ≥ 0",
            "ion_in_site": "",
            "max_coordinating_contacts": "integer ≥ 0",
            "n_ions_coordinated": "integer ≥ 0",
            "well_placed": "",
            "n_diphosphate_atoms": "integer ≥ 0",
            "diphosphate_to_cage_dist": "≥ 0 (Å)",
            "n_residues": "integer ≥ 1",
        },
    },
    "substrate_positioning": {
        "explanation": "Whether an AF3 co-folded prenyl-PP substrate is poised for catalysis, with its diphosphate at the metal cage and reactive carbon near the machinery.",
        "columns": {
            "metal_point_found": "",
            "substrate_present": "",
            "substrate_resname": "",
            "n_substrate_atoms": "integer ≥ 0",
            "substrate_plddt": "0–100",
            "diphosphate_to_cage_dist": "≥ 0 (Å)",
            "min_diphosphate_to_cage_oxygen": "≥ 0 (Å)",
            "diphosphate_to_nearest_ion": "≥ 0 (Å)",
            "reactive_carbon_to_cage_dist": "≥ 0 (Å)",
            "substrate_in_site": "",
            "n_residues": "integer ≥ 1",
        },
    },
    "motif_search": {
        "explanation": "Presence per sequence of the class-I TPS metal-binding motifs (the DDXXD aspartate-rich family and the NSE/DTE second triad).",
        "columns": {
            "DD..D": "",
            "D[DE]..[DE]": "",
            "[DE][DE]..[DE]": "",
            "(N|D)D(L|I|V).(S|T)...E": "",
        },
    },
    "soluprot": {
        "explanation": "SoluProt-predicted E. coli solubility/expressibility of each sequence; higher means more likely to express solubly.",
        "columns": {
            "soluble": "0–1",
        },
    },
    "enzyme_explorer_sequence_only": {
        "explanation": "Sequence-only EnzymeExplorer TPS classification; per-class raw scores and calibrated probabilities that a sequence is a terpene synthase.",
        "columns": {
            "2xFPP_score": "unbounded", "2xGGPP_score": "unbounded", "CPP_score": "unbounded",
            "EDSQ_score": "unbounded", "FPP_score": "unbounded", "GFPP_score": "unbounded",
            "GGPP_score": "unbounded", "GPP_score": "unbounded", "IDS_score": "unbounded",
            "TPS_score": "unbounded",
            "2xFPP_p_calibrated": "0–1", "2xGGPP_p_calibrated": "0–1", "CPP_p_calibrated": "0–1",
            "EDSQ_p_calibrated": "0–1", "FPP_p_calibrated": "0–1", "GFPP_p_calibrated": "0–1",
            "GGPP_p_calibrated": "0–1", "GPP_p_calibrated": "0–1", "IDS_p_calibrated": "0–1",
            "TPS_p_calibrated": "0–1",
        },
    },
    # --------------- COMPARATIVE / DESIGN-ONLY METRICS ---------------
    "max_sequence_identity": {
        "explanation": "Global maximum sequence identity/similarity of each design to a reference set; a novelty / near-duplicate measure.",
        "columns": {
            "sequence_identity": "0–1", "sequence_identity_hit": "",
            "sequence_similarity": "0–1", "sequence_similarity_hit": "",
        },
    },
    "min_embedding_distance": {
        "explanation": "Minimum distance in ESM-1b embedding space from each design to a reference set; an embedding-space novelty measure.",
        "columns": {
            "min_embedding_distance": "≥ 0", "min_embedding_distance_hit": "",
        },
    },
    "local_sequence_search": {
        "explanation": "Fast local (BLAST-style) best-hit sequence identity, similarity, and coverage to a reference set, plus nearest neighbours.",
        "columns": {
            "local_sequence_identity": "0–100", "local_sequence_similarity": "0–100", "local_coverage": "0–100",
        },
    },
    "swissprot_search": {
        "explanation": "DIAMOND best-hit search of each design against Swiss-Prot; flags function drift by labeling top hits as TPS vs non-TPS.",
        "columns": {
            "swissprot_top_hit": "", "swissprot_top_pident": "0–100", "swissprot_top_bitscore": "≥ 0",
            "swissprot_top_is_tps": "", "swissprot_best_nontps_pident": "0–100", "swissprot_n_tps_in_topN": "integer ≥ 0",
        },
    },
    "foldseek_swissprot_search": {
        "explanation": "Foldseek structural best-hit search of each design against AlphaFold-Swiss-Prot; the structural off-target / function-drift check.",
        "columns": {
            "foldseek_sprot_top_hit": "", "foldseek_sprot_top_tmscore": "0–1", "foldseek_sprot_top_is_tps": "",
            "foldseek_sprot_best_nontps_tmscore": "0–1", "foldseek_sprot_n_tps_in_topN": "integer ≥ 0",
        },
    },
    "structural_identity": {
        "explanation": "Foldseek best TM-score / lDDT of each design to the nearest known-TPS reference structure; structural analog of sequence identity.",
        "columns": {
            "structural_tmscore_to_known": "0–1", "structural_tmscore_to_known_hit": "",
            "structural_lddt_to_known": "0–1", "structural_lddt_to_known_hit": "",
        },
    },
    "domain_structural_identity": {
        "explanation": "Best foldseek TM-score / lDDT of each design's detected TPS domains to known reference domains, plus how many domains were found.",
        "columns": {
            "domain_structural_tmscore_to_known": "0–1", "domain_structural_tmscore_to_known_hit": "",
            "domain_structural_tmscore_to_known_type": "", "domain_structural_lddt_to_known": "0–1",
            "n_detected_domains": "integer ≥ 0",
            "domain_structural_tmscore_to_known_alpha": "0–1", "domain_structural_tmscore_to_known_beta": "0–1",
            "domain_structural_tmscore_to_known_gamma": "0–1", "domain_structural_tmscore_to_known_ids": "0–1",
            "domain_structural_tmscore_to_known_delta": "0–1", "domain_structural_tmscore_to_known_epsilon": "0–1",
            "domain_structural_tmscore_to_known_zeta": "0–1",
        },
    },
    "sdr_divergence": {
        "explanation": "Flags designs globally similar to a known TPS but divergent at specificity-determining active-site residues (the single-residue-switch regime).",
        "columns": {
            "nearest_neighbour_id": "", "nearest_neighbour_similarity": "0–1", "n_sdr_positions": "integer ≥ 0",
            "sdr_identity": "0–1", "n_sdr_mismatches": "integer ≥ 0", "specificity_divergence": "", "divergent_positions": "",
        },
    },
    "knn_label_transfer": {
        "explanation": "Predicts a coarse class for each design by a distance-weighted vote of nearest known-TPS neighbours, ensembled across sequence/embedding/structure spaces.",
        "columns": {
            "predicted_label": "", "confidence": "0–1",
        },
    },
    "substrate_class": {
        "explanation": "Predicts each design's prenyl-diphosphate substrate class by fusing a k-NN vote, the pocket-volume size band, and EnzymeExplorer scores.",
        "columns": {
            "predicted_substrate": "", "confidence": "0–1", "knn_substrate": "", "knn_confidence": "0–1",
            "pocket_volume_band": "", "substrate_agreement": "", "ee_substrate": "", "ee_score": "unbounded",
            "ee_agreement": "", "n_signals_agree": "integer 0–3", "predicted_substrate_source": "",
        },
    },
    "self_consistency": {
        "explanation": "Designability via self-consistency: whether a ProteinMPNN sequence refolds (ESMFold) back to the design; low scRMSD (<~2 Å) means designable.",
        "columns": {
            "sc_rmsd_min": "≥ 0 (Å)", "sc_rmsd_mean": "≥ 0 (Å)", "n_samples": "integer ≥ 0",
        },
    },
}
