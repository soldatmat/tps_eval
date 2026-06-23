# tps_eval tools reference

End-user documentation for every evaluation tool in `tps_eval`. Each tool follows
the repo's 3-layer shape (cluster-agnostic `scripts/run_<tool>.sh` → per-cluster
`scripts/<cluster>/jobs/<tool>.sh` → `src/<subdir>/run_<tool>.py`) and writes a
**CSV keyed by `ID`** so metrics compose/merge for filtration.

- **Sequence-branch** tools key off a FASTA and write `<input>_<tool>.csv` next to it.
- **Structure-branch** tools key off a structures **directory** and write
  `<structs_dir>_<tool>.csv` as a *sibling* of the directory. `ID` = structure
  filename stem (for AF3 `af_output` trees, the per-job name).

Tools emit **raw numbers only**; "natural TPS" bands are computed separately by the
[reference-stats pipeline](#compute_reference_stats). Conda env names and external
install paths come from [`paths.sh`](../paths.sh). The pipeline is normally run via
the [orchestrator](#run_eval_pipeline) and tool selection via
[`pipeline_tools.json`](../scripts/pipeline_tools.json).

Most tools run in the `tps_eval` conda env. Exceptions: ESMFold / ProteinMPNN /
self-consistency use `esmfold`; Aggrescan3D uses `aggrescan3d` (**Python 2.7**);
pocket descriptors use `pocket` (fpocket + P2Rank/openjdk); SoluProt uses `soluprot`;
EnzymeExplorer uses `enzyme_explorer`.

---

## Sequence tools

### motif_search
- **Purpose** — Flag presence of the two class-I TPS metal-binding motifs (the aspartate-rich DDXXD family and the NSE/DTE second triad) per sequence.
- **Inputs** — FASTA.
- **Output** — `<fasta>_motifs.csv`. Columns: `ID`, `sequence`, and one **boolean column per motif regex**. Default regexes (`scripts/run_motif_search.sh`): `DD..D` (strict), `D[DE]..[DE]`, `[DE][DE]..[DE]` (graded relaxations of DDXXD that tolerate conservative D→E substitutions), and `(N|D)D(L|I|V).(S|T)...E` (NSE/DTE). Extra motifs can be passed as positional args.
- **Method** — Compiles each motif string as a regex and records `bool(regex.search(seq))` per sequence. The relaxed DDXXD variants are nested supersets and kept as separate columns to grade hits.
- **External dependency** — none (stdlib `re` + pandas).
- **Env + source** — `tps_eval`; [`src/sequence_metrics/motif_search.py`](../src/sequence_metrics/motif_search.py). Motif definitions are centralized in [`motif_localization.py`](../src/sequence_metrics/motif_localization.py).

### motif_pair_distance
- **Purpose** — Sequence (residue) separation between the two metal-binding motifs, a coarse check that the active-site spacing is plausible.
- **Inputs** — FASTA.
- **Output** — `<fasta>_motif_pair_distance.csv`, one row per sequence. Key columns: `motif_start_distance` (= start_NSE − start_DDXXD, 1-based; positive in the canonical order), `motif_gap` (= start_NSE − end_DDXXD, the inter-motif residue gap). Helper columns record each motif's matched substring and 1-based start. Distances are NaN when either motif is absent.
- **Method** — Locates the DDXXD-family and NSE/DTE motifs via the shared `motif_localization` helper (first left-to-right match each) and reports the signed residue separations.
- **External dependency** — none.
- **Env + source** — `tps_eval`; [`src/sequence_metrics/motif_pair_distance.py`](../src/sequence_metrics/motif_pair_distance.py).

### esm_embedding
- **Purpose** — Produce per-sequence ESM-1b embeddings; an input/feeder for the embedding-distance metrics, not a standalone score.
- **Inputs** — FASTA.
- **Output** — `<fasta>_embedding_esm1b.csv`. Columns: `id` (FASTA record id) + 1280 numeric embedding-dimension columns (`0`…`1279`, the mean-pooled layer-33 representation). Note: this CSV uses `id`; the downstream distance tools re-key to `ID`.
- **Method** — Runs ESM-1b (`esm1b_t33_650M_UR50S`) and takes the mean of the final-layer (33) per-residue representations as a fixed 1280-d vector per sequence.
- **External dependency** — [ESM / ESM-1b](https://github.com/facebookresearch/esm) (Rives et al. 2021, *PNAS*). Adapted from the upstream `extract.py`.
- **Env + source** — `tps_eval`; [`src/esm/extract_embeddings.py`](../src/esm/extract_embeddings.py).

### esm_pseudo_perplexity
- **Purpose** — Sequence "naturalness" — how in-distribution a sequence is under ESM's masked language model. Lower perplexity = more natural.
- **Inputs** — FASTA.
- **Output** — `<fasta>_esm_pseudo_perplexity.csv`, keyed by `ID`. Columns: `esm_mean_pll` (mean per-residue pseudo-log-likelihood, ≤0; higher = more natural) and `esm_pseudo_perplexity` (= exp(−esm_mean_pll), ≥1; **lower = more natural**).
- **Method** — Two estimators of the per-residue pseudo-log-likelihood: `swoop` (default, fast — single unmasked forward pass, "One Fell Swoop" approximation) and `masked` (exact masked-marginal, O(L) forwards). Uses the same ESM-1b model as `esm_embedding` so naturalness is consistent with the embedding metrics.
- **External dependency** — [ESM / ESM-1b](https://github.com/facebookresearch/esm); PLL approximation following Salazar et al. 2020 (masked-LM scoring).
- **Env + source** — `tps_eval`; [`src/sequence_metrics/esm_pseudo_perplexity.py`](../src/sequence_metrics/esm_pseudo_perplexity.py).

### max_sequence_identity
- **Purpose** — Per-query maximum sequence identity to a reference set. **Self mode** (`maxid_self`) flags near-duplicates within a dataset; **gen-vs-train mode** (`maxid_gen_vs_train`) measures novelty of generated designs relative to the training set.
- **Inputs** — Query FASTA; optional reference FASTA (`--train_path`). With `--train` (self mode), each sequence is compared to all others (self-comparison excluded).
- **Output** — `<fasta>_max_sequence_identity.csv` (gen-vs-train) or `<fasta>_max_sequence_identity_self.csv` (self), keyed by `ID`. Columns: `sequence_identity`, `sequence_identity_hit` (best-matching reference id), `sequence_similarity`, `sequence_similarity_hit`. With `--top_k`, also writes a tidy `..._topk.csv` (`query_id,rank,neighbour_id,score`).
- **Method** — Global Needleman–Wunsch alignment (Biopython `PairwiseAligner`, BLOSUM62, gap open −11 / extend −1); identity and similarity are computed from the alignment and the best (max) over the reference set is reported.
- **External dependency** — [Biopython](https://biopython.org/) pairwise aligner; BLOSUM62.
- **Env + source** — `tps_eval`; [`src/sequence_metrics/max_sequence_identity.py`](../src/sequence_metrics/max_sequence_identity.py).

### min_embedding_distance
- **Purpose** — Per-query minimum ESM-embedding (Euclidean) distance to a reference set. **Self mode** (`mindist_self`) finds the nearest neighbour within a dataset; **gen-vs-train mode** (`mindist_gen_vs_train`) is an embedding-space novelty measure. Depends on `esm_embedding`.
- **Inputs** — Query embeddings CSV (`esm_embedding` output); optional reference embeddings CSV. `--train` selects self mode.
- **Output** — `<input>_embedding_esm1b_min_embedding_distance.csv` (gen-vs-train) or `..._min_embedding_distance_self.csv` (self), keyed by `ID`. Columns: `min_embedding_distance`, `min_embedding_distance_hit` (nearest reference id). With `--top_k`, also writes `..._topk.csv` (`query_id,rank,neighbour_id,score`).
- **Method** — Loads the 1280-d ESM-1b vectors, computes the pairwise Euclidean distance matrix to the reference set, and reports the minimum (and its argmin id) per query.
- **External dependency** — none beyond NumPy (consumes ESM-1b embeddings).
- **Env + source** — `tps_eval`; [`src/sequence_metrics/min_embedding_distance.py`](../src/sequence_metrics/min_embedding_distance.py).

### soluprot
- **Purpose** — Predicted *E. coli* expressibility/solubility of each sequence (sequence-based; orthogonal to the structure-based `aggregation`).
- **Inputs** — FASTA.
- **Output** — `<fasta>_soluprot.csv`, keyed by id, with SoluProt's predicted solubility score (the `soluble` target consumed by the plots). Exact column names come from the external SoluProt tool, not this repo.
- **Method** — Shells out to the external SoluProt predictor (gradient-boosted model over sequence/HMM features); needs USEARCH + TMHMM helpers and a per-job tmp dir.
- **External dependency** — [SoluProt](https://loschmidt.chemi.muni.cz/soluprot/) (Hon et al. 2021, *Bioinformatics*). Standalone install via `scripts/setup_soluprot.sh`; path/env set in `paths.sh` (`SOLUPROT_PATH`, `SOLUPROT_ENV`).
- **Env + source** — `soluprot`; wrapper [`scripts/run_soluprot.sh`](../scripts/run_soluprot.sh) (external tool, no `src/` module).

### enzyme_explorer_sequence_only
- **Purpose** — Sequence-only TPS classification — per-class probabilities that a sequence is a terpene synthase, without needing a structure.
- **Inputs** — FASTA.
- **Output** — `<fasta>_enzyme_explorer_sequence_only.csv`. Schema (from EnzymeExplorer's `predict_sequences_only` console script): `id`, `sequence`, `<class>_score`, `<class>_p_calibrated`. The plots consume the calibrated TPS probability as the `isTPS_seq` target.
- **Method** — Runs EnzymeExplorer's protein-language-model classifier (`predict_sequences_only`) with its bundled checkpoints + calibration.
- **External dependency** — [EnzymeExplorer](https://github.com/SamusRam/EnzymeExplorer) (revision branch). Installed via its own `scripts/setup_env.sh`; path/env in `paths.sh`.
- **Env + source** — `enzyme_explorer`; wrapper [`scripts/run_enzyme_explorer_sequence_only.sh`](../scripts/run_enzyme_explorer_sequence_only.sh).

### swissprot_search
- **Purpose** — Broad off-target / specificity check: what *else* does each generated sequence look like? Each top hit is labelled TPS vs non-TPS to flag function drift. Gen-only (real train TPS are trivially TPS hits).
- **Inputs** — FASTA + a DIAMOND DB built from `uniprot_sprot.fasta` (`SWISSPROT_DIAMOND_DB`) + the committed TPS-accession list.
- **Output** — `<fasta>_swissprot_search.csv`, keyed by `ID`. Columns: `swissprot_top_hit`, `swissprot_top_pident`, `swissprot_top_bitscore`, `swissprot_top_is_tps`, `swissprot_best_nontps_pident`, `swissprot_n_tps_in_topN`.
- **Method** — DIAMOND `blastp` (default `--very-sensitive`, top-N hits) vs Swiss-Prot; each hit's UniProt accession is classified TPS/non-TPS by membership in the committed TPS accession set.
- **External dependency** — [DIAMOND](https://github.com/bbuchfink/diamond) (Buchfink et al. 2021, *Nat. Methods*); [Swiss-Prot/UniProtKB](https://www.uniprot.org/). TPS set ([`src/homology_search/tps_uniprot_accessions.txt`](../src/homology_search/tps_uniprot_accessions.txt)) from the UniProt query `(reviewed:true) AND ((ec:4.2.3.*) OR (ec:5.5.1.*))`.
- **Env + source** — `tps_eval`; [`src/homology_search/swissprot_search.py`](../src/homology_search/swissprot_search.py).

### local_sequence_search
- **Purpose** — Fast **local** (BLAST-style) best-hit sequence identity/similarity to a reference set, plus top-k nearest neighbours. Complements the global full-length `max_sequence_identity` (novelty) and supplies the fast sequence-space neighbours for the k-NN/SDR tools (the Biopython all-vs-all was too slow — MMseqs2 does the 1195-seq all-vs-all in ~5 s).
- **Inputs** — query FASTA; optional reference (`--train_path`; self mode otherwise, excluding the query from its own best hit). `--backend {mmseqs2,diamond}` (default `mmseqs2`), `--top_k N`.
- **Output** — `<input>_local_sequence_search.csv`: `ID`, `local_sequence_identity` (best-hit %, both backends), `local_sequence_similarity` (DIAMOND `ppos`; NaN for mmseqs2 — `easy-search` exposes no positives field), `local_coverage`. With `--top_k`: `<input>_local_sequence_search_topk.csv` (`query_id,rank,neighbour_id,score`; score = identity %).
- **Method** — Build a DB from the reference; MMseqs2 `easy-search` (or DIAMOND `blastp`); best hit by bitscore.
- **External dependency** — [MMseqs2](https://github.com/soedinglab/MMseqs2) (bioconda) and/or [DIAMOND](https://github.com/bbuchfink/diamond); installed via `conda -c conda-forge -c bioconda --override-channels` per Aurum admin policy.
- **Env + source** — `tps_eval`; [`src/sequence_metrics/local_sequence_search.py`](../src/sequence_metrics/local_sequence_search.py).

---

## Structure tools

All structure tools accept either an AlphaFold3 `af_output` tree (reads the top-ranked `<job>/<job>_model.cif`) or a flat dir of `.pdb`/`.cif`, auto-detected via the canonical loader in [`plddt.py`](../src/structure_metrics/plddt.py). `ID` = filename stem (= AF3 job name).

### plddt
- **Purpose** — Per-structure AlphaFold/ESMFold folding-confidence summary.
- **Inputs** — Structures dir.
- **Output** — `<structs_dir>_plddt.csv`. Columns: `ID`, `mean_plddt`, `median_plddt`, `min_plddt`, `frac_plddt_confident` (fraction ≥ 70, configurable via `--confident_threshold`), `n_residues`.
- **Method** — Reads per-residue pLDDT from the **B-factor field** (PDB cols 61–66 / mmCIF `_atom_site.B_iso_or_equiv`) over protein residues (HETATM skipped) and summarizes. Valid for *predicted* structures only (for experimental structures the B-factor is a temperature factor). ESMFold writes pLDDT on a 0–1 scale that `esmfold.py` rescales ×100 so this tool reads it on the 0–100 AF convention.
- **External dependency** — Biopython parsing; consumes [AlphaFold3](https://github.com/google-deepmind/alphafold3) / [ESMFold](https://github.com/facebookresearch/esm) pLDDT.
- **Env + source** — `tps_eval`; [`src/structure_metrics/plddt.py`](../src/structure_metrics/plddt.py).

### motif_structural_distance
- **Purpose** — 3D (Å) span between the two metal-binding motifs — the structural analog of `motif_pair_distance`, approximating the active-site metal-cluster span.
- **Inputs** — Structures dir.
- **Output** — `<structs_dir>_motif_structural_distance.csv`. Columns: `ID`, `motif_centroid_distance` (Å, between the centroids of the DDXXD- and NSE/DTE-coordinating Cα atoms — primary metric), `motif_min_ca_distance` (closest Cα–Cα approach), `n_residues`. NaN when a motif is absent.
- **Method** — Derives the 1-letter sequence from the structure, locates both motifs with the shared `motif_localization` helper, and measures distances between their coordinating-residue Cα atoms.
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/structure_metrics/motif_structural_distance.py`](../src/structure_metrics/motif_structural_distance.py).

### active_site_geometry
- **Purpose** — Side-chain-level geometry of the catalytic carboxylate cage — whether the metal-coordinating oxygens actually converge on a single locus that could hold the trinuclear Mg²⁺/Mn²⁺ cluster. Apo-robust (no metals/ligand modelled needed).
- **Inputs** — Structures dir; optional `--templates ID[,ID...]` (reference IDs present in the dir, e.g. `1ps1`,`5eat`) to enable the constellation-RMSD columns.
- **Output** — `<structs_dir>_active_site_geometry.csv`. Columns: `ID`, `carboxylate_convergence_radius` (RMS distance of coordinating O atoms from their centroid; ~6–9 Å for a competent site), `n_coordinating_oxygens`, `metal_point_void` (clearance Å at the oxygen centroid; a clash/filled centroid is bad), `catalytic_constellation_rmsd` (best superposition RMSD to a reference constellation; only with `--templates`), `best_template`, `n_residues`. Geometry is NaN when a motif or coordinating oxygen is absent.
- **Method** — Locates the motifs, gathers the carboxylate/hydroxyl side-chain oxygens (Asp OD1/OD2, Glu OE1/OE2, Asn OD1, Ser OG, Thr OG1), and computes convergence radius + centroid clearance. Optional Biopython `Superimposer` RMSD of the matched Cα+Cβ constellation to reference templates.
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/structure_metrics/active_site_geometry.py`](../src/structure_metrics/active_site_geometry.py).

### radius_of_gyration
- **Purpose** — Global shape/compactness descriptors of the fold.
- **Inputs** — Structures dir.
- **Output** — `<structs_dir>_radius_of_gyration.csv`. Columns: `ID`, `radius_of_gyration` (Å, unweighted over Cα), `asphericity` (Å²; 0 = spherical, large = elongated), `acylindricity` (Å²), `principal_radius_1/2/3` (Å, √λ of the gyration-tensor eigenvalues), `n_residues`. Raw geometric numbers, no expected band.
- **Method** — Computes Rg = √(mean‖rᵢ−r_com‖²) over the Cα atoms and the gyration-tensor eigen-decomposition (λ1≥λ2≥λ3) for shape anisotropy.
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/structure_metrics/radius_of_gyration.py`](../src/structure_metrics/radius_of_gyration.py).

### pocket_descriptors
- **Purpose** — Characterize the *catalytic* cavity (the one holding the metal cluster + prenyl-diphosphate substrate); catalytic-pocket volume is the headline "molecular-ruler ↔ product chain-length" signal.
- **Inputs** — Structures dir; optional `--fpocket` / `--prank` executables (P2Rank via `P2RANK_PATH`; omit to skip the P2Rank cross-check).
- **Output** — `<structs_dir>_pocket_descriptors.csv`. Columns: `ID`, `metal_point_found`, fpocket: `catalytic_pocket_volume` (Å³), `pocket_hydrophobicity`, `pocket_enclosure`, `pocket_n_alpha_spheres`, `pocket_total_sasa`, `pocket_depth`, `pocket_sasa_per_volume`, `fpocket_catalytic_pocket_found`; P2Rank: `p2rank_catalytic_site_score`, `p2rank_catalytic_pocket_rank`, `p2rank_catalytic_pocket_found`; `n_residues`. Descriptors are NaN when no detected pocket coincides with the metal point (itself a red flag) or a motif is absent. `pocket_sasa_per_volume` (DERIVED) = `pocket_total_sasa / catalytic_pocket_volume` (Å⁻¹), a specific-surface-area shape/compactness descriptor; note it is size-dependent (~1/radius) so it tracks cavity size rather than isolating shape (a size-free sphericity would need a matched cavity surface+volume, which fpocket's lining-atom SASA and alpha-sphere volume are not).
- **Method** — Anchors on the carboxylate-cage metal point (reuses `active_site_geometry` + the shared motif localizer), then selects the fpocket pocket (Voronoi alpha-spheres) and the P2Rank pocket nearest/enclosing that point and reports each engine's descriptors.
- **Reproducibility** — `catalytic_pocket_volume` is fpocket's Monte-Carlo volume estimate and is **stochastic ~1.5% run-to-run** (median rel. diff 1.5%, p95 4%, max ~8% across the 1348-protein MARTS-DB set) — fpocket does not fix the MC seed. `pocket_total_sasa`, `pocket_depth`, alpha-sphere counts, hydrophobicity/enclosure, and the P2Rank scores are deterministic and reproduce exactly. So `pocket_sasa_per_volume` inherits the same ~1.5% volume noise in its denominator. Treat the volume band as ±a few percent; don't expect bit-identical volumes when re-running.
- **External dependency** — [fpocket](https://github.com/Discngine/fpocket) (Le Guilloux et al. 2009, *BMC Bioinformatics*); [P2Rank](https://github.com/rdk/p2rank) (Krivák & Hoksza 2018, *J. Cheminform.*; prebuilt release at `P2RANK_PATH`).
- **Env + source** — `pocket` (fpocket + openjdk for P2Rank's `prank`); [`src/structure_metrics/pocket_descriptors.py`](../src/structure_metrics/pocket_descriptors.py).

### domain_composition
- **Purpose** — TPS structural-domain composition per design: how many and which TPS domain types are present.
- **Inputs** — Structures dir (`.pdb` files); or `--detections_json` to parse an existing EE sidecar instead of re-detecting.
- **Output** — `<structs_dir>_domain_composition.csv`, one row per design (zero-domain designs included). Columns: `ID`, `n_domains`, per-type counts `n_alpha`, `n_beta`, `n_gamma`, `n_ids`, `n_delta`, `n_epsilon`, `n_zeta`, and `domain_architecture` (hyphen-joined type string, `""` for zero-domain). Consumed by `plot_domains`.
- **Method** — EnzymeExplorer's `detect_domains` aligns each structure (PyMOL + foldseek; CPU-only, no PLM) against seven curated TPS domain templates. Designs with zero detected domains are enumerated independently from the dir and left-joined so every design gets exactly one row.
- **External dependency** — [EnzymeExplorer](https://github.com/SamusRam/EnzymeExplorer) (domain detector); foldseek; PyMOL.
- **Env + source** — `enzyme_explorer`; [`src/enzyme_explorer/domain_composition.py`](../src/enzyme_explorer/domain_composition.py).

### aggregation
- **Purpose** — Structure-based aggregation propensity (an expressibility signal orthogonal to the sequence-based `soluprot`).
- **Inputs** — Structures dir.
- **Output** — `<structs_dir>_aggregation.csv`, keyed by `ID`. Columns: `a3d_avg_score`, `a3d_total_score`, `a3d_max_score`, `a3d_min_score`, `a3d_total_pos_score`, `n_residues`. Positive = aggregation-prone surface-exposed hydrophobic; per-structure A3D failures emit a NaN row.
- **Method** — Runs Aggrescan3D in **static mode** (never the slow dynamic CABS-flex mode), which scores spatially-clustered surface hydrophobic patches per residue, then reduces the per-residue A3D output to per-ID scalars.
- **External dependency** — [Aggrescan3D / A3D](https://github.com/lcbio/aggrescan3d) (Kuriata et al. 2019, *NAR*), vendored at `vendor/aggrescan3d` (**Python 2.7**).
- **Env + source** — `aggrescan3d` (Py2.7 — module is Py2-compatible); [`src/structure_metrics/aggregation.py`](../src/structure_metrics/aggregation.py).

### foldseek_swissprot_search
- **Purpose** — Broad *structural* off-target check: each design's nearest AlphaFold-Swiss-Prot structures, top hit labelled TPS vs non-TPS (the structural analog of `swissprot_search`).
- **Inputs** — Structures dir + the foldseek AlphaFold/Swiss-Prot DB (`AFDB_SWISSPROT_DB`) + the committed TPS-accession list.
- **Output** — `<structs_dir>_foldseek_swissprot_search.csv`, keyed by `ID`. Columns: `foldseek_sprot_top_hit`, `foldseek_sprot_top_tmscore`, `foldseek_sprot_top_is_tps`, `foldseek_sprot_best_nontps_tmscore`, `foldseek_sprot_n_tps_in_topN`.
- **Method** — foldseek `easy-search` (top-N, TM-score), each AFDB target accession (`AF-<ACC>-F1-model_v4` → `<ACC>`) classified TPS/non-TPS by the committed accession set.
- **External dependency** — [foldseek](https://github.com/steineggerlab/foldseek) (van Kempen et al. 2024, *Nat. Biotechnol.*); [AlphaFold-Swiss-Prot DB](https://alphafold.ebi.ac.uk/).
- **Env + source** — `tps_eval`; [`src/homology_search/foldseek_swissprot_search.py`](../src/homology_search/foldseek_swissprot_search.py).

### structural_identity
- **Purpose** — Foldseek structural identity (best TM-score / lddt) of each design to the nearest *known TPS* reference structure — the structural analog of `max_sequence_identity`. Requires a reference-structures dir (`--known_structs_dir`).
- **Inputs** — Generated structures dir + known-TPS reference structures dir.
- **Output** — `<structs_dir>_structural_identity.csv`, keyed by `ID`. Columns: `structural_tmscore_to_known`, `structural_tmscore_to_known_hit`, `structural_lddt_to_known`, `structural_lddt_to_known_hit`. With `--top_k`, also writes `..._structural_identity_topk.csv` (`query_id,rank,neighbour_id,score`).
- **Method** — foldseek all-vs-known structural alignment; per query, keeps the best (max) `alntmscore` and `lddt` and the matching reference stem (self-hits excluded for self-search).
- **External dependency** — [foldseek](https://github.com/steineggerlab/foldseek).
- **Env + source** — `tps_eval`; [`src/structure_metrics/run_structural_identity.py`](../src/structure_metrics/run_structural_identity.py) (alignment in [`src/foldseek/structure_alignment.py`](../src/foldseek/structure_alignment.py)).

### domain_structural_identity
- **Purpose** — Structural identity at the **domain** level rather than the whole chain: detect each design's TPS structural domains (α/β/γ/δ/ε/ζ/IDS) and foldseek-align *each domain* to the curated known martsDB reference domains. Catches designs whose overall fold drifts but whose individual catalytic/support domains still match a known TPS (and vice-versa), and reports how many domains were detected.
- **Inputs** — Generated structures dir. Reference DOMAIN structures default to EnzymeExplorer's curated set (`$ENZYME_EXPLORER_PATH/data/detected_domains/martsDB_detected_domains/domains`); override with `--known_domain_structures_root`.
- **Output** — `<structs_dir>_domain_structural_identity.csv`, keyed by `ID`. Columns: `domain_structural_tmscore_to_known`, `domain_structural_tmscore_to_known_hit`, `domain_structural_tmscore_to_known_type`, `domain_structural_lddt_to_known`, `n_detected_domains`, and per-type bests `domain_structural_tmscore_to_known_{alpha,beta,gamma,ids,delta,epsilon,zeta}`.
- **Method** — EnzymeExplorer's `detect_domains` carves each design into its constituent domains, then foldseek aligns the detected domains against the reference domains; per design, keeps the best TM-score/lddt overall and per reference domain-type.
- **External dependency** — [EnzymeExplorer](https://github.com/) `detect_domains` + [foldseek](https://github.com/steineggerlab/foldseek) (both live in the `enzyme_explorer_prod` env on Aurum).
- **Env + source** — `enzyme_explorer_prod`; [`src/structure_metrics/run_domain_structural_identity.py`](../src/structure_metrics/run_domain_structural_identity.py) (reuses [`src/foldseek/domain_alignment.py`](../src/foldseek/domain_alignment.py)).

### proteinmpnn_score
- **Purpose** — Sequence-given-fold likelihood: how compatible each design's *own* sequence is with its backbone. Lower = more likely given the fold.
- **Inputs** — Structures dir.
- **Output** — `<structs_dir>_proteinmpnn_score.csv`. Columns: `ID`, `proteinmpnn_nll` (mean per-residue NLL over all residues = ProteinMPNN `global_score`), `proteinmpnn_score_designed` (NLL over designed residues = `score`). For a fully-designed monomer these coincide.
- **Method** — Shells out to the vendored ProteinMPNN once per structure with `--score_only 1` (scores the native sequence read from the PDB), reading back the mean NLL from the `score_only/*.npz` arrays.
- **External dependency** — [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) (Dauparas et al. 2022, *Science*), vendored at `vendor/ProteinMPNN` (ships its own weights).
- **Env + source** — `esmfold` (`PROTEINMPNN_ENV` defaults to it); [`src/structure_metrics/proteinmpnn_score.py`](../src/structure_metrics/proteinmpnn_score.py).

### self_consistency
- **Purpose** — Designability via self-consistency scRMSD: does there exist a sequence that ProteinMPNN likes for this fold *and* that ESMFold refolds back to the same shape? **Heavy / opt-in** (default off; GPU, N folds per structure).
- **Inputs** — Structures dir; `--num_seqs` (default 8), `--ids`/`--limit` to validate on a few structures.
- **Output** — `<structs_dir>_self_consistency.csv`. Columns: `ID`, `sc_rmsd_min` (best of N), `sc_rmsd_mean`, `n_samples` (folds that succeeded). A design is self-consistent/designable when `sc_rmsd_min` < ~2 Å.
- **Method** — Per backbone: sample N sequences with ProteinMPNN → refold each with ESMFold → Cα-align each refold to the original design (Biopython `Superimposer`) and record RMSD.
- **External dependency** — [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) + [ESMFold](https://github.com/facebookresearch/esm).
- **Env + source** — `esmfold` (has torch + transformers + Biopython; ProteinMPNN shelled out with the same python); [`src/structure_metrics/self_consistency.py`](../src/structure_metrics/self_consistency.py).

### aromatic_lining
- **Purpose** — Aromatic (Trp/Tyr/Phe) residues lining the catalytic pocket that stabilize carbocation intermediates — a cyclization-capability proxy. Counts are apo-robust; ring orientation is softer on open apo sites.
- **Inputs** — structures dir; `--cutoff` (default ~10 Å around the metal point).
- **Output** — `<structs_dir>_aromatic_lining.csv`: `ID`, `n_pocket_aromatics`, `n_trp`/`n_tyr`/`n_phe`, `aromatic_fraction`, `n_inward_facing_aromatics`, `metal_point_found`. NaN/0 when the metal point can't be located.
- **Method** — Locates the active-site metal point (the relaxed `active_site_geometry.metal_point`), takes pocket residues within the cutoff, counts aromatics and the rings whose face points into the cavity within cation-π range.
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/structure_metrics/aromatic_lining.py`](../src/structure_metrics/aromatic_lining.py).

### diphosphate_sensor
- **Purpose** — Whether the design supplies the basic residues (Arg/Lys) and the conserved RY pair that anchor/ionize the substrate diphosphate at the metal site. Presence is apo-robust; exact rotamers are softer.
- **Inputs** — structures dir; `--cutoff` (~12 Å), `--ry_dist` (~6 Å).
- **Output** — `<structs_dir>_diphosphate_sensor.csv`: `ID`, `metal_point_found`, `n_diphosphate_basic_residues`, `n_arg`, `n_lys`, `has_RY_pair`, `n_RY_pairs`, `n_residues`.
- **Method** — Relaxed metal point; counts Arg/Lys whose terminal N atoms point toward the metal/diphosphate region within the cutoff; detects an Arg+Tyr (RY) pair near the site.
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/structure_metrics/diphosphate_sensor.py`](../src/structure_metrics/diphosphate_sensor.py).

### ion_site_check
- **Purpose** — Geometric validation that the Mg²⁺/Mn²⁺ ions **AlphaFold3 co-folded** (any holo `--af3_cofold mg*`) actually land in the catalytic carboxylate cage. AF3 free-ion placement is a *hypothesis*; every other active-site tool is apo-robust and anchors on the protein-derived cage point only. This is the one tool that READS the ion HETATMs (which the shared parser otherwise skips) and checks them against that expected cage. **Only carries signal for AF3 holo folds** — apo / ESMFold structures have no ions and report a graceful not-applicable row (`n_ions_modelled=0`, distances NaN, bools False).
- **Inputs** — structures dir; `--site_radius` (default 5.0 Å, in-site distance from the cage centroid), `--coord_cutoff` (default 2.8 Å, Mg–O coordination; real bonds ~2.0–2.5 Å), `--min_coord_contacts` (default 2), `--ion_resnames` (default `MG MN`), `--diphosphate_resnames` (default `POP PPV PPK`).
- **Output** — `<structs_dir>_ion_site_check.csv`: `ID`, `metal_point_found`, `n_ions_modelled`, `min_ion_to_cage_dist`, `n_ions_in_site`, `ion_in_site`, `max_coordinating_contacts`, `n_ions_coordinated`, `well_placed`, `mg_canonical_motif_coordination`, `n_motif_coord_asp`, `n_motif_coord_nse`, `mg_to_motif_dist`, `n_diphosphate_atoms`, `diphosphate_to_cage_dist`, `n_residues`. Distance columns are NaN and bool columns False when there are no ions or no metal point.
- **Method** — Computes the expected apo cage point via the canonical relaxed `active_site_geometry.metal_point` (centroid of the DDXXD + NSE/DTE coordinating side-chain oxygens), reads the ion HETATMs, and measures: nearest ion → cage distance, ions inside the site sphere, and per-ion coordinating-oxygen contacts at Mg–O bonding distance (`well_placed` = ≥1 ion coordinated by ≥`min_coord_contacts` cage oxygens — the strict validation). For `mg_ppi`, also the diphosphate-centroid → cage distance. **Reference-independent addition:** `mg_canonical_motif_coordination` flags whether the modelled ions are coordinated by the ACTUAL DDXXD + NSE/DTE motif carboxylates (scanning all motif occurrences via the shared `sequence_metrics.motif_localization`), which — unlike the apo `metal_point` — does not mislocalize for two-domain folds; `n_motif_coord_asp`/`n_motif_coord_nse` count the coordinating residues per motif and `mg_to_motif_dist` is the ion-centroid → DDXXD-oxygen-centroid distance.
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/structure_metrics/ion_site_check.py`](../src/structure_metrics/ion_site_check.py). For the trinuclear Mg²⁺ cluster geometry see Christianson, D. W. *Chem. Rev.* **2017**, *117*, 11570–11648.

### substrate_positioning
- **Purpose** — Is the **AlphaFold3 co-folded prenyl-PP substrate** (`--af3_cofold mg_<sub>` forced, or `mg_ee` per-design) actually POISED FOR CATALYSIS — its diphosphate at the DDXXD/NSE metal cage and its reactive carbon (C1, where the carbocation forms once PPi leaves) held near the catalytic machinery rather than in bulk solvent? Like `ion_site_check`, it READS the ligand HETATMs. **Only carries signal for SUBSTRATE holo folds** — apo / ESMFold / Mg-only / `mg_ppi` (bare pyrophosphate, no substrate) report a graceful not-applicable row (`substrate_present=False`, geometry NaN).
- **Inputs** — structures dir; `--site_radius` (default 6.0 Å, informational diphosphate-centroid → cage band), `--coord_cutoff` (default 4.0 Å, the robust in-cage test: closest diphosphate atom → cage carboxylate oxygen), `--ion_resnames` (default `MG MN`), `--min_substrate_carbons` (default 5; a HETATM residue with ≥1 P and ≥ this many C is the prenyl-PP substrate — distinguishes it from `POP` and the ions), `--substrate_resname` (force a specific ligand residue name instead of auto-detect).
- **Output** — `<structs_dir>_substrate_positioning.csv`: `ID`, `metal_point_found`, `substrate_present`, `substrate_resname`, `n_substrate_atoms`, `substrate_plddt` (mean ligand B-factor/pLDDT), `diphosphate_to_cage_dist` (centroid→centroid; inherits the shared metal_point, can be inflated by a splayed cage — see note), `min_diphosphate_to_cage_oxygen`, `diphosphate_to_nearest_ion`, `diphosphate_to_ion_centroid`, `reactive_carbon_to_cage_dist`, `reactive_carbon_to_nearest_ion`, `reactive_carbon_to_ion_centroid`, `substrate_in_site`, `n_residues`.
- **Method** — Auto-detects the substrate ligand by composition (most carbons among residues with ≥1 P and ≥`min_substrate_carbons` C); its diphosphate = the ligand's P + O atoms, its reactive carbon = the ligand carbon nearest the diphosphate. Reuses the canonical `active_site_geometry.metal_point` + coordinating-oxygen set. `substrate_in_site` is based on the **oxygen** distance (`min_diphosphate_to_cage_oxygen` ≤ `coord_cutoff`), NOT the centroid distance, because the relaxed coordinating-oxygen centroid can sit far from the diphosphate even when it is well coordinated. **Reference-independent addition:** `diphosphate_to_ion_centroid` / `reactive_carbon_to_nearest_ion` / `reactive_carbon_to_ion_centroid` measure the same poising against the COFOLDED Mg cluster instead of the apo `metal_point` (which mislocalizes for two-domain folds). **Caveat:** AF3 ligand placement from SMILES is a *hypothesis* — this tool is exactly the downstream check.
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/structure_metrics/substrate_positioning.py`](../src/structure_metrics/substrate_positioning.py). Substrate SMILES in [`src/alphafold/cofold_substrates.py`](../src/alphafold/cofold_substrates.py). Christianson, D. W. *Chem. Rev.* **2017**, *117*, 11570–11648.

### cyclization_geometry
- **Purpose** — Cheap, **reference-independent**, *necessary-not-sufficient* geometric signals that the AF3 co-folded prenyl-PP substrate is organized to CYCLIZE (beyond merely being bound, which `substrate_positioning` covers): (1) **substrate fold** — is the prenyl chain curled so a distal carbon can reach C1 (first ring-closure geometry) vs splayed out? and (2) **cation-π track** — are aromatic side chains lined along the substrate carbons to stabilize the migrating carbocation cascade? **Only carries signal for SUBSTRATE holo folds** — apo / ESMFold / Mg-only / `mg_ppi` report a graceful not-applicable row (`substrate_present=False`, geometry NaN).
- **Inputs** — structures dir; `--aromatic_cutoff` (default 6.0 Å, substrate-carbon → aromatic-ring-centroid cation-π contact), `--farchain_bonds` (default 6, bonds from C1 to count a carbon as distal for fold-back), `--ion_resnames` (default `MG MN`), `--min_substrate_carbons` (default 5), `--substrate_resname` (force a ligand resname).
- **Output** — `<structs_dir>_cyclization_geometry.csv`: `ID`, `substrate_present`, `n_substrate_carbons`, `substrate_rgyr` (ligand-heavy-atom radius of gyration; smaller = more compact/folded), `foldback_c1_to_distal` (min distance C1 → a chain carbon ≥`farchain_bonds` bonds away; small = curled toward C1 = cyclization-compatible), `substrate_endtoend`, `n_aromatic_carbon_contacts`, `frac_aromatic_track`, `n_aromatics_lining`, `mean_carbon_to_aromatic`, `n_residues`.
- **Method** — Reuses `substrate_positioning.read_substrate_ligand` (validated AF3/Boltz2 ligand+ion parse). C1 = the ligand carbon nearest the diphosphate; chain topology is a BFS over C–C bonds (<1.8 Å). Aromatic ring centroids are taken from PHE/TYR/TRP/HIS side chains. All metrics are intrinsic to the ligand or substrate-carbon → protein-aromatic, so they never use the apo `metal_point`. **Caveat:** NECESSARY-NOT-SUFFICIENT — a folded, aromatic-lined substrate is consistent with a competent cyclase but does NOT establish the product (the carbocation cascade, hydride/methyl shifts and quench are unmodeled).
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/structure_metrics/cyclization_geometry.py`](../src/structure_metrics/cyclization_geometry.py). Christianson, D. W. *Chem. Rev.* **2017**, *117*, 11570–11648.

### global_confidence
- **Purpose** — Whole-fold confidence (pTM / iPTM), complementing per-residue pLDDT.
- **Inputs** — `--pae_dir` (per-structure `<ID>_pae.npz` carrying `ptm`/`iptm`, saved at fold time by ESMFold or the AF3 extract_pae step); optional `--structs_dir` for output naming.
- **Output** — `<structs_dir>_global_confidence.csv`: `ID`, `ptm`, `iptm` (when present). NaN when the npz/ptm is missing.
- **Method** — Reads the `ptm`/`iptm` scalars from the saved PAE npz.
- **External dependency** — numpy.
- **Env + source** — `tps_eval`; [`src/structure_metrics/global_confidence.py`](../src/structure_metrics/global_confidence.py). PAE/pTM retention lives in [`src/esmfold/esmfold.py`](../src/esmfold/esmfold.py) + [`src/alphafold/extract_pae.py`](../src/alphafold/extract_pae.py).

### interdomain_pae
- **Purpose** — Confidence in the **relative orientation** of the design's domains — a multi-domain failure mode pLDDT can't see.
- **Inputs** — structures dir + `--pae_dir` (`<ID>_pae.npz`).
- **Output** — `<structs_dir>_interdomain_pae.csv`: `ID`, `mean_interdomain_pae`, `max_interdomain_pae`, `n_domains` (+ optional per-domain-pair columns). N/A for single-domain designs or missing PAE.
- **Method** — EnzymeExplorer `detect_domains` gives per-domain residue ranges; reduces the off-diagonal inter-domain PAE blocks (averaged both directions, since PAE is asymmetric).
- **External dependency** — EnzymeExplorer (domain detector) + numpy.
- **Env + source** — `enzyme_explorer_prod`; [`src/structure_metrics/interdomain_pae.py`](../src/structure_metrics/interdomain_pae.py).

---

## Folding (structure producers)

These produce the `structs/` dir of `<ID>.pdb` consumed unchanged by the structure tools. They are **not** wired into the orchestrator (v2); run them first, then pass `--structs_dir`.

### alphafold3
- **Purpose** — Fold (and optionally co-fold with ligands/ions) sequences with AlphaFold3. **Aurum-only.** **Orchestrator-wired**: pass `--fold alphafold3` to `run_eval_pipeline.py` and it fans the gen FASTA out into one AF3 job per sequence, extracts PAE, then runs the whole structure branch on the result (no pre-supplied `--structs_dir` needed).
- **Inputs** — Standalone: a CSV of proteins (+ optional ligand SMILES / ion CCD codes); see `--protein_id_column_names` / `--protein_sequence_column_names` / `--ligand_*` / `--ion_*`. Via the orchestrator: just the gen FASTA — the fan-out wrapper (`scripts/run_alphafold_fanout.sh`) converts it to a CSV. **Co-fold** the class-I TPS active site for a HOLO prediction with `--af3_cofold`: `none` (default, apo protein only); `mg` (the trinuclear Mg²⁺ cluster — 3× CCD `MG`, ligated by DDXXD + NSE/DTE); `mg_ppi` (the cluster + a bare diphosphate head group, CCD `POP`/pyrophosphate²⁻, a substrate-agnostic stand-in); `mg_gpp` / `mg_fpp` / `mg_ggpp` / `mg_gfpp` (the cluster + ONE forced prenyl-PP substrate as SMILES — `cofold_substrates.py` — for EVERY design); or `mg_ee` (the cluster + each design's *own* EnzymeExplorer-predicted substrate; the fan-out splits designs into per-substrate groups + a Mg-only fallback for non-co-foldable EE calls). With `--fold alphafold3` the pipeline **auto-chains** EE → cofold: it submits the sequence branch, then a small `eval_pipeline_continuation` job (`afterok` on it) that re-invokes the pipeline once the EE CSV exists — so one command runs EE *and* substrate co-folding (the login-node fold driver can't itself wait on the `ee_seq` job, so the lightweight continuation does). Pass `--enzymeexplorer_csv` only when folding is not in-pipeline (pre-computed EE) or to override. Any non-`none` mode enables the holo tools (`ion_site_check`, `substrate_positioning`); `none` / `--no_holo_tools` turns co-folding and those tools off. Output filenames stay `<ID>.pdb` regardless. **Caveat:** AF3 free-ion/ligand placement is a *hypothesis*, not ground truth — verify the Mg/diphosphate land at the DDXXD/NSE cage (`ion_site_check` / `substrate_positioning`) before trusting the holo geometry (Christianson 2017; AF3 co-folding literature).
- **Output** — Under the orchestrator, a `<gen>_af3/` work dir holding `af_output/` (per-job AF3 trees), `structs/<ID>.pdb` (CIF→PDB extracted, pLDDT in the B-factor via the patched `vendor/cif_to_pdb`), and `pae/<ID>_pae.npz` (from the `extract_pae` step). The structure tools auto-detect the layout.
- **Method** — A login-node driver (`run_alphafold_jobs.py`) submits one AF3 SLURM job per sequence (custom `b32_128_gpu --constraint=alphafold3` partition), skipping existing structures, and prints the N job ids; the orchestrator captures them so the structure branch `afterok`-waits on all N. Each job folds + extracts CIF→PDB; a following `extract_pae` job (`src/alphafold/extract_pae.py`) populates the PAE dir. PAE-consumers (`global_confidence`, `interdomain_pae`) wait on that extraction step.
- **External dependency** — [AlphaFold3](https://github.com/google-deepmind/alphafold3) (Abramson et al. 2024, *Nature*).
- **Env + source** — `tps_eval` (for the driver); [`src/alphafold/run_alphafold_jobs.py`](../src/alphafold/run_alphafold_jobs.py), fan-out wrapper [`scripts/run_alphafold_fanout.sh`](../scripts/run_alphafold_fanout.sh), PAE extraction [`src/alphafold/extract_pae.py`](../src/alphafold/extract_pae.py). See README "Running AlphaFold".

### esmfold
- **Purpose** — Fold single-chain sequences with ESMFold — a fast single-sequence alternative to AF3 that runs on **both clusters**. **Orchestrator-wired**: pass `--fold esmfold` to `run_eval_pipeline.py` and it folds the generated FASTA first, then runs the whole structure branch on the result (no pre-supplied `--structs_dir` needed).
- **Inputs** — FASTA; `--save_dir` (structs out dir), `--pae_dir` (PAE out dir; default `<save_dir>_pae/`), `--chunk_size`/`--device` tuning, `--no-skip_existing`, `--no-save_pae`.
- **Output** — One `<ID>.pdb` per FASTA record in the structs dir (ID = record id = filename stem), mirroring the AlphaFold `structs/` layout, plus one `<ID>_pae.npz` per record in the PAE dir (consumed by `global_confidence` + `interdomain_pae`). Per-residue pLDDT is written to the B-factor field, rescaled 0–1 → 0–100 so `plddt` reads it.
- **Method** — Runs `facebook/esmfold_v1` (HuggingFace transformers); sequences > ~600 aa trigger chunked attention to bound GPU memory. When driven by the orchestrator, the pipeline derives `<gen>_esmfold_structs/` + `<gen>_esmfold_structs_pae/` and makes every structure Step depend on the `esmfold_gen` producer.
- **External dependency** — [ESMFold](https://github.com/facebookresearch/esm) / `facebook/esmfold_v1` (Lin et al. 2023, *Science*).
- **Env + source** — `esmfold`; [`src/esmfold/esmfold.py`](../src/esmfold/esmfold.py) (wrapper `scripts/run_esmfold.sh`).

---

## Representation / embedding producers

These emit a **per-`id` feature CSV** (first column `id`, then feature dims), keyed by the structure/sequence stem (`Enzyme_marts_ID` for the MARTS-DB set). They are **standalone — NOT orchestrator `Step`s** (not wired into `run_eval_pipeline.py` / `pipeline_tools.json`); run them directly, then feed the CSV to [`run_visualization`](#run_visualization) (or any downstream consumer). They produce *representations* of a whole protein set, not per-design eval metrics.

### saprot_embedding
- **Purpose** — SaProt-650M **structure-aware** per-protein embeddings: a PLM representation that fuses sequence and structure (each residue's amino acid letter is concatenated with its foldseek 3Di structural-state token), mean-pooled to one vector per structure. A structure-conditioned alternative to the sequence-only `esm_embedding`.
- **Inputs** — `--structs_dir` (dir of `<ID>.pdb`; stem == row key), `--foldseek` (binary path, used to derive the 3Di tokens), `--output_csv`; optional `--ids_csv`/`--id_column` (restrict to a reference id set + report coverage), `--saprot_repo` (clone providing `utils.foldseek_util.get_struc_seq`), `--model_location`, `--chain` (default `A`), `--truncation_seq_length` (default 1024), `--nogpu`.
- **Output** — Feature CSV keyed by `id` (PDB filename stem): first column `id`, then embedding dims `0..1279` (mean-pooled layer-output over real residues, excluding BOS/EOS). Failed structures are skipped (reported, not NaN rows).
- **Method** — For each PDB, SaProt's `get_struc_seq` shells out to foldseek to build the structure-aware (SA) token sequence (AA + 3Di interleaved; `plddt_mask="auto"` masks low-pLDDT residues only for self-identifying AlphaFold PDBs, so ESMFold PDBs are unmasked), then runs `EsmForMaskedLM` and mean-pools `last_hidden_state` over real residues. Mirrors `src/esm/extract_embeddings.py` conventions (first column `id`, then dims `0..D-1`).
- **External dependency** — [SaProt](https://github.com/westlake-repl/SaProt) (Su et al., *ICLR* 2024), model `westlake-repl/SaProt_650M_AF2` (ESM-2 650M backbone, vocab = AA × 3Di); [foldseek](https://github.com/steineggerlab/foldseek) (3Di tokens); HuggingFace transformers.
- **Env + source** — `saprot` (separate env; defaults point at the Karolina shared-project install); [`src/saprot/extract_saprot_embeddings.py`](../src/saprot/extract_saprot_embeddings.py) (wrapper [`scripts/run_saprot_embedding.sh`](../scripts/run_saprot_embedding.sh)).

### ee_domain_features
- **Purpose** — EnzymeExplorer **domain-comparison** features: the structure/function feature block of EE's production model (`PlmDomainsRandomForest__tps_esm-1v-subseq_..._domains_subset`). Each detected TPS structural domain is compared to EE's curated reference functional-domain modules by foldseek TM-score, yielding a per-protein `1 − TM-score` profile over the production `domains_subset` modules.
- **Inputs** — `--sequences_csv` (defines the output row set via `--id_column`, default `Enzyme_marts_ID`), `--structs_dir` (dir of `<id>.pdb`), `--output_csv`; optional `--scratch_dir`, `--n_jobs` (default 16), `--classifier_pkl` (production fold-classifier bundle defining the reference-module columns), `--reuse_cached` (skip detection+comparison, rebuild the matrix from cached pickles).
- **Output** — Feature CSV keyed by `id`: first column `id`, then one column per reference module (`<known_module_id>`) — value = `1 − max TM-score` over the protein's detected domains of the matching type against that module (`1.0` where no comparison exists). The column universe is the UNION of selected modules across the 5 production fold-classifiers (≈897 cols = the `domains_subset` selection). First/second detected α-domain are split into `alpha1`/`alpha2` exactly as EE's classifier does.
- **Method** — Detect TPS structural domains per structure (EE `domain_detections`), compare each to the curated reference modules with foldseek (`comparing_to_known_domains_foldseek`, restricted to `domains_subset.pkl`), then build `1 − dom_feat` per `easy_predict.py`. Detection/comparison run on **sanitized underscore-free IDs** (e.g. `marts_E00000` → `martsE00000`) to dodge EE's `query.split('_')[0]` keying bug, then map back to the original `id`. CPU-only.
- **External dependency** — [EnzymeExplorer](https://github.com/SamusRam/EnzymeExplorer) (Samusevich et al., *bioRxiv* 2024.01.29.577750) `domain_detections` + `comparing_to_known_domains_foldseek`; foldseek; PyMOL.
- **Env + source** — `enzyme_explorer` (run from the EE `scripts/` dir so `data/` resolves to the production bundle); [`src/enzyme_explorer/extract_ee_domain_features.py`](../src/enzyme_explorer/extract_ee_domain_features.py).

### ee_esm1v_embeddings
- **Purpose** — EnzymeExplorer **ESM-1v-TPS** PLM embeddings: the PLM feature block of EE's production model (the companion of `ee_domain_features`). The exact mean-pooled representation EE's `easy_predict.py` feeds the classifier, from the TPS-finetuned "subseq" ESM-1v checkpoint.
- **Inputs** — `--sequences_csv` (with `--id_column`, default `Enzyme_marts_ID`, and `--sequence_column`, default `Aminoacid_sequence`), `--output_csv`; optional `--batch_size` (default 8), `--max_seq_len` (default 1022).
- **Output** — Feature CSV keyed by `id`: first column `id`, then `emb_0..emb_1279` (1280-d, mean-pooled layer-33 over residue tokens, excluding BOS/EOS). The production classifier's `plm_feat_indices_subset` is the full `range(1280)`, so this raw embedding IS the PLM block — no subsetting.
- **Method** — Loads ESM-1v (`esm1v_t33_650M_UR90S_1`) with EE's TPS-finetuned subseq checkpoint via `get_model_and_tokenizer("esm-1v-finetuned-subseq")`, runs `compute_embeddings` (layer 33, mean-pooled, `max_len`-truncated exactly as `easy_predict.py`'s domain branch). Sequences are read from the input CSV (covers proteins without a structure). GPU strongly recommended.
- **External dependency** — [EnzymeExplorer](https://github.com/SamusRam/EnzymeExplorer) (Samusevich et al., *bioRxiv* 2024.01.29.577750); [ESM-1v](https://github.com/facebookresearch/esm); torch.
- **Env + source** — `enzyme_explorer` (run from the EE `scripts/` dir so `data/plm_checkpoints/...` resolves); [`src/enzyme_explorer/extract_ee_esm1v_embeddings.py`](../src/enzyme_explorer/extract_ee_esm1v_embeddings.py).

### active_site_features
- **Purpose** — A 32-d **active-site / cation-specific-residue property + geometry profile** of the catalytic pocket of each class-I TPS, intended for a class-coloured landscape map of first-cyclization specificity. Premise (Durairaj et al., *PLOS Comput. Biol.* 2021): TPS first-cyclization specificity is dominated by the active-site contour residues lining the Mg²⁺ carboxylate cage, not the global fold. Alignment-free (no multiple alignment across the 22 classes / multiple folds).
- **Inputs** — `--structs_dir` (dir of `<id>.pdb` ESMFold structures), `--marts_csv` (`TPS_first_cyclization.csv`, for the class label merge), `--output`; optional `--radius` (shell radius, default 12 Å), `--exclude_class` (default `"12"` — drops the OSC outlier **by label**; pass `""` to keep it), `--features_only_csv` (also dump the raw structure-keyed feature CSV before the class merge). The lower-level `extract_active_site_features.active_site_features_dir` can run on a structs dir alone.
- **Output** — CSV keyed by `id` (== `Enzyme_marts_ID` == structure stem). Metadata: `product_class_id`, `product_class_marts_id`, `substrate_name`, `n_product_classes`, `all_product_class_ids`, `metal_point_found`, `n_shell_residues`, `n_residues`, `radius_A`. **32 feature columns**: 7 property-group fractions (`frac_aromatic/aliphatic/acidic/basic/polar/glycine/proline`) + 20 per-AA fractions (`frac_aa_<X>`) + 5 geometry (`carboxylate_convergence_radius`, `n_coordinating_oxygens`, `shell_radius_of_gyration`, `mean_dist_to_metal_point`, `n_aromatic_within_8A`). Feature columns are **NaN** for any protein with no locatable metal point (`metal_point_found=False`). Raw numbers only — no class-conditional banding.
- **Method** — Anchors on the carboxylate-cage **metal point** (centroid of the DDXXD + NSE/DTE coordinating side-chain oxygens) via the canonical `structure_metrics.active_site_geometry.metal_point` (the same definition used by `aromatic_lining` / `pocket_descriptors`), selects the active-site **shell** (residues with any atom within `--radius` of the metal point), and featurizes it as an order-/length-invariant property + composition profile plus cage geometry. OSC class-12 enzymes (different fold, no Mg²⁺ cluster, no real DDXXD/NSE-DTE) are excluded by **label** because the relaxed motif regex can spuriously match — `metal_point_found` alone does not flag them.
- **External dependency** — Biopython (read-only reuse of `active_site_geometry` + the shared motif localizer); [MARTS-DB](https://marts-db.org/) first-cyclization labels.
- **Env + source** — `tps_eval`; [`src/specificity/extract_active_site_features.py`](../src/specificity/extract_active_site_features.py) (logic) + [`run_extract_active_site_features.py`](../src/specificity/run_extract_active_site_features.py) (argv + MARTS-DB class merge).

---

## Function (structure-dependent)

### enzyme_explorer
- **Purpose** — TPS classification with structures (per-class scores, richer than the sequence-only variant). **Not yet wired into the orchestrator (v2).**
- **Inputs** — A FASTA *or* a sequences CSV (`ID`,`sequence`) **plus** a `--structs_dir`.
- **Output** — An `<input>_enzyme_explorer/` output *directory* (the revision-branch `predict_with_structures` schema — no longer a single `_enzyme_explorer.csv`). The plots consume the TPS probability as the `isTPS` target.
- **Method** — Runs EnzymeExplorer's structure-aware predictor (`predict_with_structures`), combining the PLM classifier with structural domain features.
- **External dependency** — [EnzymeExplorer](https://github.com/SamusRam/EnzymeExplorer) (revision branch).
- **Env + source** — `enzyme_explorer`; wrapper [`scripts/run_enzyme_explorer.sh`](../scripts/run_enzyme_explorer.sh).

---

## Aggregator & visualization

### plots
- **Purpose** — The comparison aggregator: merge every enabled metric CSV (by `ID`) across datasets and render the comparison figures. Effectively always on in the orchestrator unless excluded.
- **Inputs** — One or more FASTAs with `--data_names` / `--data_colors` (the metric CSVs are discovered next to each FASTA / structs dir by naming convention); optional `--targets`, `--save_dir`.
- **Output** — Plot images written to `--save_dir` (a `plots/` dir beside the gen FASTA in the orchestrator). Numeric metrics → boxplot + density; categorical/boolean metrics (motif presence, `domain_architecture`, `*_top_is_tps`) → count plots. Targets and their CSV sources are mapped in `src/plot/constants.py`.
- **Method** — Loads each target's source CSV(s), merges per dataset, and draws per-target comparison panels.
- **External dependency** — matplotlib (+ pandas).
- **Env + source** — `tps_eval`; [`src/plot/run_plots.py`](../src/plot/run_plots.py) / [`plot_comparison.py`](../src/plot/plot_comparison.py); targets in [`src/plot/constants.py`](../src/plot/constants.py).

### plot_domains
- **Purpose** — Render per-design PyMOL images of detected TPS domains overlaid on the full structure (visual companion to `domain_composition`). Standalone (not orchestrator-wired).
- **Inputs** — A selection CSV (`--structures_column_name`, default `ID`), `--structures_root`, `--domain_structures_root`, `--domains_pkl` (EnzymeExplorer detected-domains pickle), `--output_root`.
- **Output** — Per-design PNGs (and optional `.pse` PyMOL sessions) under `--output_root`.
- **Method** — Loads each full structure and its per-domain PDBs in PyMOL, colors/overlays the detected domains, and ray-traces images.
- **External dependency** — PyMOL (`pymol-open-source` by default; see README "licensed PyMOL Incentive").
- **Env + source** — `tps_eval`; `src/pymol/plot_domains.py` (wrapper [`scripts/run_plot_domains.sh`](../scripts/run_plot_domains.sh)).

### plot_residue_similarity
- **Purpose** — For each (query, matched-known) pair, render PyMOL images coloring the new structure by per-residue sequence similarity (BLOSUM90) to its matched known structure, plus a side-by-side alignment (visual companion to `structural_identity`). Standalone.
- **Inputs** — A selection CSV (defaults: query column `query`, known column `max_alntmscore_target`), `--structures_root`, `--known_structures_root`, `--output_root`.
- **Output** — Per-pair PNGs (and optional `.pse` sessions) under `--output_root`.
- **Method** — Aligns each query to its matched known structure, scores per-residue BLOSUM90 similarity, and colors/ray-traces the structure plus an alignment view in PyMOL.
- **External dependency** — PyMOL; BLOSUM90.
- **Env + source** — `tps_eval`; `src/pymol/plot_residue_similarity.py` (wrapper [`scripts/run_plot_residue_similarity.sh`](../scripts/run_plot_residue_similarity.sh)).

### run_visualization
- **Purpose** — Dataset-level **landscape-map** visualization: project a whole protein set (MARTS-DB-style) into 2D from any representation and render it coloured by first-cyclization class, with a carbon-ordered substrate-type palette so the ESM / SaProt / EE / active-site / similarity maps are directly comparable. Standalone (not orchestrator-wired). Distinct from `plots` (per-design metric-comparison charts) and the `src/pymol/` 3D renders — this lays out a *dataset* in 2D representation space.
- **Inputs** — One input mode plus one label mode. **Feature mode** `--features CSV` (first column `id`, then numeric feature dims — e.g. a SaProt / EE / active-site producer CSV); methods `pca,tsne,umap,pacmap`. **Pairs mode** `--pairs TSV` (headerless all-vs-all similarity table; `--pair-cols q,t,...`, `--sim-col`; distance = 1 − sim, missing pairs → 1); methods `umap,tsne,pcoa` (precomputed-distance). Labels (for colour) via `--label-col COL` (label already in the features CSV), `--labels-parallel CSV --class-col COL` (labels by ROW POSITION, row-aligned with the features CSV), or `--labels-join CSV --id-col COL --class-col COL` (join by id; first class per id). `--methods` (comma list, one panel each), `--title` (required), `--output` (required), `--exclude-cols`, `--footnote`.
- **Output** — A multi-panel **PNG figure** (one scatter panel per method, 200 dpi) at `--output` — **not** a CSV. Points coloured by first-cyclization class (22-class palette); a grouped legend (by substrate type) is parked outside the data area. Rows with any NaN feature are dropped; PCA/PCoA panel titles carry the %-variance of the two axes.
- **Method** — `dimensionality_reduction.py` backends: PCA (z-scored SVD), t-SNE (PCA-init from top-50 PCs, or `metric=precomputed`), UMAP (euclidean or `metric=precomputed`), PaCMAP (feature only), PCoA / classical MDS (precomputed distance only). `landscape_map.render_panels` draws the class-coloured scatter using `palette.make_palette` (each substrate type = one hue family — mono C10 Greens → sesqui C15 Blues → di C20 Reds → sester C25 Purples → sterol/triterpene brown singleton — shaded by class within the family). t-SNE/UMAP/PaCMAP are imported lazily.
- **External dependency** — numpy, pandas, matplotlib; scikit-learn (t-SNE), [umap-learn](https://github.com/lmcinnes/umap), [pacmap](https://github.com/YingfanWang/PaCMAP) (only for the requested methods).
- **Env + source** — `tps_eval`; [`scripts/run_visualization.py`](../scripts/run_visualization.py) (runner) over [`src/visualization/`](../src/visualization/) (`dimensionality_reduction.py`, `landscape_map.py`, `palette.py`).

---

## Label transfer

### knn_label_transfer
- **Purpose** — Predict a **coarse class** for each generated design by a distance-weighted vote of its nearest MARTS-DB known-TPS neighbours, **ensembled across the three similarity spaces** (`max_sequence_identity`, `min_embedding_distance`, `structural_identity`), with an honest **leave-one-out calibration** on MARTS-DB. **Label-agnostic:** the class assignments are an INPUT (`--label_file`, a `reference_id,label` CSV) — swap the file to change the labeling (first-cyclization class, size class, substrate, …); nothing in the logic hardcodes a particular labeling.
- **Inputs** — The per-design **top-k CSVs** emitted by the three tools' `--top_k` flag (`<input>_local_sequence_search_topk.csv`, `<input>_min_embedding_distance_topk.csv`, `<structs_dir>_structural_identity_topk.csv`; each `query_id,rank,neighbour_id,score`; the pipeline feeds the sequence space from the fast `local_sequence_search`), a `--label_file`, and a `--calibration` JSON. Any space may be omitted (the design abstains there). Two subcommands: `calibrate` (consumes the MARTS-DB **self** top-k CSVs → calibration JSON) and `predict` (consumes the **design** top-k CSVs + calibration JSON → predictions).
- **Output** — `predict`: CSV keyed by `ID` with `predicted_label`, `confidence` (calibrated), and per space `predicted_label_<space>`, `conf_<space>`, `nn_similarity_<space>`. Designs below τ in **all** spaces **abstain** (`predicted_label = "unknown"`, `confidence = 0`) — novel designs should land here. `calibrate`: a committable JSON artifact (`src/reference_stats/knn_calibration_<labeling>.json`) with per-space + ensemble accuracy, the chosen τ, and the binned nn_similarity→P(correct) calibration curve.
- **Method** — Per space: convert each neighbour's score to a similarity in [0,1] (`identity%/100`; TM-score as-is; embedding distance → `1/(1+d)`), strip the foldseek `_<chain>` suffix from structural `neighbour_id` (only when the stripped stem is a known label id) before joining to the label file, **ignore neighbours below a space-specific τ** (abstain if none qualify), distance-weight the vote, and normalize to a per-class posterior (argmax = predicted). Confidence = `winning_fraction × top-k_agreement × nearest-neighbour_similarity`, reported **calibrated** via the LOO curve. Ensemble = average of the per-space posteriors (each space's argmax contributes only when it does not abstain). Calibration: leave-one-out over the labeled MARTS-DB set (self top-k **excluding self**), measuring accuracy vs nearest-neighbour similarity per space and ensembled; τ is the lowest nn_similarity at which empirical P(correct) ≥ `--target_accuracy` (default 0.5), floored at the literature prior (≈40 % identity for class transfer, TM≈0.5 fold floor; embedding has no prior so it is purely empirical).
- **External dependency** — none for the transfer itself (pandas/numpy); the top-k CSVs come from the three existing tools. The first-cyclization label file is derived from the companion `tps-first-cyclization-knn` table via [`make_first_cyclization_labels.py`](../src/knn/make_first_cyclization_labels.py).
- **Env + source** — `tps_eval`; [`src/knn/knn_label_transfer.py`](../src/knn/knn_label_transfer.py) (logic) + [`run_knn_label_transfer.py`](../src/knn/run_knn_label_transfer.py) (argv) + [`scripts/run_knn_label_transfer.sh`](../scripts/run_knn_label_transfer.sh).

### sdr_divergence
- **Purpose** — Flags designs that are **globally close to a known-product TPS but diverge at the specificity-determining active-site residues** — the TEAS/HPS single-residue-switch regime that global-similarity transfer misses. The companion negative-filter to the k-NN.
- **Inputs** — structures dir + `--known_structs_dir` + the sequence/structural `--*_topk` neighbour CSVs (the nearest known-TPS neighbour); optional `--sdr_panel <file>` of explicit specificity positions (else a structure-derived active-site panel around the metal point).
- **Output** — `<structs_dir>_sdr_divergence.csv`: `ID`, `nearest_neighbour_id`, `nearest_neighbour_similarity`, `n_sdr_positions`, `sdr_identity`, `n_sdr_mismatches`, `specificity_divergence` (bool), `divergent_positions`.
- **Method** — Take the rank-1 neighbour from the top-k, superpose the design onto it (Biopython), compare residues at the SDR/active-site positions, and flag high global similarity + low SDR-residue identity.
- **External dependency** — Biopython.
- **Env + source** — `tps_eval`; [`src/specificity/sdr_divergence.py`](../src/specificity/sdr_divergence.py).

### substrate_class
- **Purpose** — Predict each design's **prenyl-diphosphate substrate class** (GPP/C10 mono, FPP/C15 sesqui, GGPP/C20 di, …) by **fusing three independent signals**: (1) the label-agnostic k-NN vote run with a *substrate* label file + substrate calibration (the call of record), (2) the `pocket_descriptors` catalytic-pocket volume mapped to a coarse size band (the active-site "molecular ruler"), and (3) EnzymeExplorer's per-substrate sequence-only scores (argmax substrate). Reports the fused call plus how many of the corroborating signals agree.
- **Inputs** — At least one of the three `--{sequence,embedding,structural}_topk` neighbour CSVs (same feeders as the k-NN), a substrate `--label_file` (`reference_id,label`, e.g. [`src/knn/substrate_labels.csv`](../src/knn/substrate_labels.csv)) + `--calibration` (substrate calibration JSON); optional `--pocket_csv` (pocket_descriptors output) and `--enzymeexplorer_csv` (enzyme_explorer_sequence_only output) for the two cross-checks. The orchestrator wires all of these when `--train_path`, `--structs_dir`, and `--known_structs_dir` are present.
- **Output** — `<input>_substrate_class.csv`, keyed by `ID`. Columns: `predicted_substrate`, `confidence`, `knn_substrate`, `knn_confidence`, `pocket_volume_band`, `substrate_agreement`, `ee_substrate`, `ee_score`, `ee_agreement`, `n_signals_agree`, `predicted_substrate_source`.
- **Method** — Runs `knn_label_transfer.transfer_labels` over the substrate labels for the primary call, maps the fpocket volume to a size band (coarse, monotonic with chain length), takes EE's argmax substrate, and records within-one-size-class agreement of each cross-check with the k-NN call. When k-NN abstains but EE is confident, falls back to EE's argmax (flagged via `predicted_substrate_source`).
- **Env + source** — `tps_eval`; [`src/knn/run_substrate_class.py`](../src/knn/run_substrate_class.py) (logic in [`src/knn/substrate_class.py`](../src/knn/substrate_class.py)).

---

## Reference & orchestration

### run_eval_pipeline
- **Purpose** — Cluster-agnostic declarative orchestrator. Submits every enabled tool's SLURM job in dependency order, skips steps whose output already exists (idempotent/resumable), and chains deps as a single `--dependency=afterok:…`. Supersedes the per-cluster `submit_all.sh`.
- **Inputs/usage** — `python scripts/run_eval_pipeline.py --cluster <aurum|karolina> --fasta_path gen.fasta [--train_path train.fasta] [--fold esmfold|alphafold3 | --structs_dir structs/] [--known_structs_dir known/] [--self_consistency] [--dry-run]`. `--fold` produces the structures first (`esmfold` on both clusters; `alphafold3` Aurum-only, a per-sequence fan-out) so you don't need `--structs_dir`.
- **Tool selection** — Driven by [`pipeline_tools.json`](#pipeline_tools) (each key has a `default` on/off + `branch` + one-line `description`). CLI overrides, in precedence order: `--only A,B` (run only these, + plots), `--include A,B` (force-enable), `--exclude A,B` (force-disable), `--list-tools` (print the catalog and exit). `--self_consistency` is a back-compat alias for `--include self_consistency`.
- **Scope** — Full sequence branch + plots + the structure-consuming branch (everything that *reads* structures) + both structure producers: the **ESMFold producer** (`--fold esmfold`, both clusters) and the **AlphaFold3 per-sequence fan-out** (`--fold alphafold3`, Aurum-only; a login-node driver submits one AF3 job per sequence and the engine waits on all N). **Not yet ported (v2):** `enzyme_explorer`-with-structures. (AF3 holo co-folding of the Mg²⁺ cluster + diphosphate is available via `--af3_cofold`.)
- **Env + source** — pure stdlib (runs on a login node, no conda env); [`scripts/run_eval_pipeline.py`](../scripts/run_eval_pipeline.py).

### pipeline_tools
- **Purpose** — The config that drives orchestrator tool selection: maps each tool KEY → `{default, branch, description}`. The one-liners here are what `--list-tools` renders and what this doc's table descriptions mirror.
- **Notes** — If absent/stale, the orchestrator falls back to the built-in `DEFAULT_TOOLS` table (and fills in any missing keys from it). To register a new tool: add one entry here + tag its `Step(s)` with `tool="<key>"` in `build_steps`.
- **Source** — [`scripts/pipeline_tools.json`](../scripts/pipeline_tools.json).

### compute_reference_stats
- **Purpose** — The reference-stats pipeline: compute the "natural TPS" bands. Runs the *intrinsic-property* metric tools on the MARTS-DB known-TPS reference set and aggregates each metric column into summary statistics. **Standalone — deliberately NOT in `run_eval_pipeline.py`.**
- **Inputs/usage** — `scripts/compute_reference_stats.sh --cluster <c> --fasta_path data/train/TPS_sequences.fasta --ref_dir <outside-repo dir> [--structs_dir <MARTS-DB structs>] [--sequence_only] [--aggregate_only]`.
- **Output** — `src/reference_stats/marts_db_metric_stats.json` (committable). Per metric/column: numeric columns get count/mean/std/min/percentiles (p1…p99)/max; categorical/boolean columns get a frequency table. Included metrics: `motif_pair_distance`, `esm_pseudo_perplexity`, `plddt`, `motif_structural_distance`, `active_site_geometry`, `aggregation`, `domain_composition`, `proteinmpnn_score`, `radius_of_gyration`. Excluded (inherently comparative): identity / embedding / structural distance / broad searches / scRMSD.
- **Method** — Submits the existing per-tool job scripts on the staged MARTS-DB inputs, then chains [`aggregate_reference_stats.py`](../src/reference_stats/aggregate_reference_stats.py) (column-type-driven, metric-agnostic) as a dependent job.
- **External dependency** — [MARTS-DB](https://marts-db.org/) (the natural-TPS reference set).
- **Env + source** — `tps_eval` (aggregator); [`scripts/compute_reference_stats.sh`](../scripts/compute_reference_stats.sh) + [`src/reference_stats/aggregate_reference_stats.py`](../src/reference_stats/aggregate_reference_stats.py).
