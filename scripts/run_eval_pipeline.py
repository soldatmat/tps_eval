#!/usr/bin/env python3
"""Cluster-agnostic, declarative orchestrator for the tps_eval pipeline.

Replaces scripts/<cluster>/submit_all.sh. The pipeline is defined ONCE as a list
of Steps (name, job script, args, output path, dependencies, capability gate);
the engine submits each step's SLURM job in dependency order, skipping any step
whose output already exists (idempotent/resumable) and chaining dependencies as a
single `--dependency=afterok:id1:id2:...` (the old bash emitted multiple
`--dependency=` flags, of which SLURM honored only the last).

Cluster differences live in CLUSTERS (submit command, dependency/log formatting,
capabilities) — to add a cluster, add an entry, not a forked script.

Scope: covers the SEQUENCE branch + plots (motif, esm, min/max distance/identity,
soluprot, enzyme_explorer_sequence_only) AND the structure-level metrics that
consume already-folded structures: pLDDT (folding confidence) and foldseek
structural identity to the nearest known TPS (--structs_dir / --known_structs_dir).
Structure PRODUCERS are wired via --fold: 'esmfold' (both clusters, one whole-FASTA
job) and 'alphafold3' (Aurum-only, a per-sequence FAN-OUT — one AF3 job per design via
a login-node driver that captures the N job ids, then a PAE-extraction step). Either
folds the generated FASTA into a structs dir (+ PAE) the whole structure branch then
consumes, so no pre-supplied --structs_dir is needed. What is NOT yet ported (v2):
EnzymeExplorer-with-structures (AF3 holo co-folding via --af3_cofold IS wired).

Usage:
    python scripts/run_eval_pipeline.py --cluster aurum \
        --fasta_path gen.fasta [--train_path train.fasta] \
        [--fold esmfold|alphafold3 | --structs_dir structs/] [--known_structs_dir known_tps_structs/] \
        [--train_embeddings_path <csv>] [--data_colors dodgerblue goldenrod] [--dry-run]

Config-driven tool selection
-----------------------------
Which metric tools run is decided by a tools config — `scripts/pipeline_tools.json`
(loaded via stdlib `json`; if absent/unreadable, the built-in DEFAULT_TOOLS table is
used). It maps each tool KEY -> {default: bool, branch: "sequence"|"structure",
description: str}. A KEY may back several Steps (e.g. `motif` -> motif_gen +
motif_train); enabling/disabling the key affects all of them. Each Step is tagged
with its `tool` key so the engine can filter. `self_consistency` defaults OFF (heavy
scRMSD); `plots` is the aggregator and is effectively always on unless excluded.

CLI overrides (applied over the config defaults, in this precedence order):
    --only A,B,...      run ONLY these tool keys (+ plots unless also excluded).
    --include A,B,...   force-enable on top of defaults (e.g. --include self_consistency).
    --exclude A,B,...   force-disable.
    --list-tools        print every tool key (default/branch/description) and exit.
Unknown tool keys are rejected with the valid list. `--self_consistency` is kept as a
back-compat alias for `--include self_consistency`.

To register a NEW tool (e.g. pocket_descriptors): add one entry to
scripts/pipeline_tools.json (and mirror it in DEFAULT_TOOLS below for the
config-missing fallback), then add its `out_<tool>()` helper + a Step tagged
`tool="<key>"` in build_steps. The engine/CLI pick it up with no further changes.

Pure stdlib — runs on a login node, no conda env needed (it only shells out to
`sbatch` and checks output paths). DO NOT add a non-stdlib dependency (no pyyaml);
the config is JSON on purpose.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # scripts/.. == repo root
TOOLS_CONFIG = os.path.join(REPO, "scripts", "pipeline_tools.json")

# Single pipeline-wide top-k for the three feeder tools (local_sequence_search,
# min_embedding_distance, structural_identity). Their *_topk.csv outputs feed the
# k-NN label-transfer and SDR-divergence consumers. Built-in fallback; overridden by
# the JSON config "_settings.top_k" and then by the --top_k CLI flag.
DEFAULT_TOP_K = 15

# Default k-NN / substrate reference artifacts (committed under src/, see CLAUDE.md
# "/data/ is gitignored" -> committable reference artifacts live under src/).
#
# REFERENCE = MARTS-DB 2026-06-12, STRUCTURE SOURCE = ESMFold. The structural-space
# calibration consumes structural top-k neighbours, so a default must commit to one
# fold source; we default to ESMFold (full PAE coverage, one consistent fold source on
# both clusters). The AF3-space variants live alongside as knn_calibration_*_af3.json --
# point the calibration args at them for an AF3 reference (note: ion_site_check only
# carries signal on AF3 holo folds, not ESMFold/apo). Labels are structure-independent
# (no _esmfold/_af3 split). Release + structure are also recorded inside the band JSONs
# (marts_db_{esmfold,af3}_metric_stats.json -> structure_source / marts_db_release).
DEFAULT_KNN_LABEL_FILE = os.path.join(REPO, "src", "knn", "first_cyclization_labels.csv")
DEFAULT_KNN_CALIBRATION = os.path.join(
    REPO, "src", "reference_stats", "knn_calibration_first_cyclization_esmfold.json")

# Default substrate-class combiner reference artifacts (the substrate label file is the
# MARTS `Type`-derived GPP/FPP/GGPP/... map; its own calibration mirrors the k-NN one but
# over the substrate label vocabulary).
DEFAULT_SUBSTRATE_LABEL_FILE = os.path.join(REPO, "src", "knn", "substrate_labels.csv")
DEFAULT_SUBSTRATE_CALIBRATION = os.path.join(
    REPO, "src", "reference_stats", "knn_calibration_substrate_esmfold.json")


# --------------------------------------------------------------------------- #
# Cluster adapters                                                            #
# --------------------------------------------------------------------------- #
def _slurm_dep(job_ids: List[str]) -> List[str]:
    # SLURM ANDs colon-separated job ids in a SINGLE --dependency flag.
    return [f"--dependency=afterok:{':'.join(job_ids)}"] if job_ids else []


def _slurm_log(log_dir: str) -> List[str]:
    return [f"--output={log_dir}/%x.%j.out"]


def _slurm_jobid(submit_stdout: str) -> str:
    m = re.search(r"Submitted batch job (\d+)", submit_stdout)
    if not m:
        raise RuntimeError(f"could not parse a SLURM job id from:\n{submit_stdout}")
    return m.group(1)


def _slurm_fanout_ids(driver_stdout: str) -> List[str]:
    """Parse the job ids a fan-out driver submitted. run_alphafold_jobs.py prints
    ``AlphaFold job IDs: ['123', '456']`` as its LAST line (its own contract). Extract
    every integer on that line. Empty (e.g. all structures already present) -> [], so
    downstream steps simply don't wait on a fold job."""
    m = re.search(r"AlphaFold job IDs:\s*(.*)", driver_stdout)
    if not m:
        raise RuntimeError(
            f"could not find the 'AlphaFold job IDs:' line in fan-out output:\n{driver_stdout}")
    return re.findall(r"\d+", m.group(1))


CLUSTERS: Dict[str, dict] = {
    "aurum": {"submit": "sbatch", "dep": _slurm_dep, "log": _slurm_log,
              "jobid": _slurm_jobid, "fanout_ids": _slurm_fanout_ids, "caps": {"alphafold"}},
    "karolina": {"submit": "sbatch", "dep": _slurm_dep, "log": _slurm_log,
                 "jobid": _slurm_jobid, "fanout_ids": _slurm_fanout_ids, "caps": set()},  # AF3 Aurum-only
}


# --------------------------------------------------------------------------- #
# Output-path conventions (must match the tools' own _get_save_path naming)   #
# --------------------------------------------------------------------------- #
def resolve_sbatch_account(cluster: str, explicit: Optional[str] = None) -> str:
    """Resolve the SLURM account and export it as $SBATCH_ACCOUNT so every submitted
    job inherits it (this orchestrator builds `sbatch` directly and runs it via
    subprocess, which inherits os.environ -- sbatch reads $SBATCH_ACCOUNT natively).

    Precedence: --account > an SBATCH_ACCOUNT already in the environment >
    scripts/<cluster>/config.sh (which sources the uncommitted, per-install
    config.local.sh where the project grant id lives -- kept out of git because it is
    install-specific and changes on allocation renewal). Returns "" if none resolved;
    the caller warns, since a cluster's default account often lacks a partition
    association and would get jobs rejected.
    """
    if explicit:
        os.environ["SBATCH_ACCOUNT"] = explicit
        return explicit
    acct = os.environ.get("SBATCH_ACCOUNT", "").strip()
    cfg = os.path.join(REPO, "scripts", cluster, "config.sh")
    if not acct and os.path.isfile(cfg):
        try:
            acct = subprocess.run(
                ["bash", "-c", f'. "{cfg}" >/dev/null 2>&1; printf %s "${{SBATCH_ACCOUNT:-}}"'],
                capture_output=True, text=True, check=False,
            ).stdout.strip()
        except Exception:  # noqa: BLE001 - best-effort; caller warns if still empty
            acct = ""
    if acct:
        os.environ["SBATCH_ACCOUNT"] = acct
    return acct


def _base(fasta: str) -> str:
    for ext in (".fasta", ".fa"):
        if fasta.endswith(ext):
            return fasta[: -len(ext)]
    return os.path.splitext(fasta)[0]


def out_motif(f): return _base(f) + "_motifs.csv"
def out_esm(f): return _base(f) + "_embedding_esm1b.csv"
def out_mindist(f): return _base(f) + "_embedding_esm1b_min_embedding_distance.csv"
def out_mindist_self(f): return _base(f) + "_embedding_esm1b_min_embedding_distance_self.csv"
def out_maxid(f): return _base(f) + "_max_sequence_identity.csv"
def out_maxid_self(f): return _base(f) + "_max_sequence_identity_self.csv"
# top-k feeder outputs (consumed by knn / sdr_divergence). The feeder tools key
# the top-k CSV off the SAME input as their single-best CSV:
#   max_sequence_identity  -> <fasta>_max_sequence_identity_topk.csv
#   min_embedding_distance -> <embeddings_csv>_min_embedding_distance_topk.csv
#   structural_identity    -> <structs_dir>_structural_identity_topk.csv
def out_maxid_topk(f): return _base(f) + "_max_sequence_identity_topk.csv"
def out_mindist_topk(emb_csv): return os.path.splitext(emb_csv)[0] + "_min_embedding_distance_topk.csv"
# local_sequence_search (MMseqs2/DIAMOND): the FAST local counterpart of
# max_sequence_identity. Keyed off the fasta. Its _topk.csv is the SEQUENCE-space
# feeder for knn + sdr_divergence (replacing the slow Biopython max_sequence_identity
# for that role; max_sequence_identity stays as the plain GLOBAL novelty metric).
def out_local_search(f): return _base(f) + "_local_sequence_search.csv"
def out_local_search_topk(f): return _base(f) + "_local_sequence_search_topk.csv"
# Self mode (within-set novelty) writes to a distinct _self path (passed explicitly via
# --save_path) so it does not collide with the gen-vs-train output, which shares the
# default <input>_local_sequence_search.csv name. Mirrors maxid_self vs maxid.
def out_local_search_self(f): return _base(f) + "_local_sequence_search_self.csv"
def out_soluprot(f): return _base(f) + "_soluprot.csv"
def out_ee_seq(f): return _base(f) + "_enzyme_explorer_sequence_only.csv"
def out_motif_pair(f): return _base(f) + "_motif_pair_distance.csv"
def out_swissprot_search(f): return _base(f) + "_swissprot_search.csv"
def out_esm_ppl(f): return _base(f) + "_esm_pseudo_perplexity.csv"
# Structure-branch outputs are keyed by the structures DIRECTORY, not the fasta:
# the tools save "<structs_dir>_<tool>.csv" as a sibling of the directory.
def out_plddt(d): return d.rstrip(os.sep) + "_plddt.csv"
def out_structural_identity(d): return d.rstrip(os.sep) + "_structural_identity.csv"
def out_structural_identity_topk(d): return d.rstrip(os.sep) + "_structural_identity_topk.csv"
# Consumer (dependency-heavy) outputs keyed by the gen fasta / structs dir:
def out_knn(f): return _base(f) + "_knn_label_transfer.csv"
def out_sdr_divergence(d): return d.rstrip(os.sep) + "_sdr_divergence.csv"
def out_motif_struct(d): return d.rstrip(os.sep) + "_motif_structural_distance.csv"
def out_active_site_geom(d): return d.rstrip(os.sep) + "_active_site_geometry.csv"
def out_radius_of_gyration(d): return d.rstrip(os.sep) + "_radius_of_gyration.csv"
def out_pocket_descriptors(d): return d.rstrip(os.sep) + "_pocket_descriptors.csv"
def out_domain_composition(d): return d.rstrip(os.sep) + "_domain_composition.csv"
def out_aggregation(d): return d.rstrip(os.sep) + "_aggregation.csv"
def out_foldseek_swissprot(d): return d.rstrip(os.sep) + "_foldseek_swissprot_search.csv"
def out_proteinmpnn(d): return d.rstrip(os.sep) + "_proteinmpnn_score.csv"
def out_self_consistency(d): return d.rstrip(os.sep) + "_self_consistency.csv"
def out_interdomain_pae(d): return d.rstrip(os.sep) + "_interdomain_pae.csv"
def out_global_confidence(d): return d.rstrip(os.sep) + "_global_confidence.csv"
def out_aromatic_lining(d): return d.rstrip(os.sep) + "_aromatic_lining.csv"
def out_diphosphate_sensor(d): return d.rstrip(os.sep) + "_diphosphate_sensor.csv"
def out_ion_site_check(d): return d.rstrip(os.sep) + "_ion_site_check.csv"
def out_substrate_positioning(d): return d.rstrip(os.sep) + "_substrate_positioning.csv"
def out_domain_structural_identity(d): return d.rstrip(os.sep) + "_domain_structural_identity.csv"
# substrate_class is keyed off the gen FASTA (it fuses sequence + structure signals).
def out_substrate_class(f): return _base(f) + "_substrate_class.csv"


# --------------------------------------------------------------------------- #
# Tool catalog (config-driven on/off)                                         #
# --------------------------------------------------------------------------- #
# Built-in fallback used when scripts/pipeline_tools.json is missing/unreadable.
# Keep in sync with that file. To register a NEW tool, add ONE entry here (and to
# the JSON), then tag its Step(s) with tool="<key>" in build_steps -- e.g.:
#   "pocket_descriptors": {"default": True, "branch": "structure",
#                          "description": "Pocket descriptors ..."},
DEFAULT_TOOLS: Dict[str, dict] = {
    "motif":                {"default": True,  "branch": "sequence",  "description": "DDXXD / NSE-DTE motif presence search."},
    "motif_pair":           {"default": True,  "branch": "sequence",  "description": "Sequence distance between the two metal-binding motifs."},
    "esm":                  {"default": True,  "branch": "sequence",  "description": "ESM-1b embeddings (feeds mindist_*)."},
    "esm_ppl":              {"default": True,  "branch": "sequence",  "description": "ESM pseudo-perplexity (sequence likelihood)."},
    "maxid_self":           {"default": True,  "branch": "sequence",  "description": "Max pairwise sequence identity within the dataset."},
    "local_sequence_search": {"default": True, "branch": "sequence",  "description": "Fast LOCAL (MMseqs2) sequence identity/similarity/coverage search; self (within-set novelty) + gen-vs-train. Its _topk.csv is the k-NN/SDR sequence-space feeder."},
    "mindist_self":         {"default": True,  "branch": "sequence",  "description": "Min ESM-embedding distance within the dataset (needs esm)."},
    "soluprot":             {"default": True,  "branch": "sequence",  "description": "SoluProt predicted solubility."},
    "ee_seq":               {"default": True,  "branch": "sequence",  "description": "EnzymeExplorer sequence-only TPS classification."},
    "swissprot_search":     {"default": True,  "branch": "sequence",  "description": "DIAMOND search vs Swiss-Prot (gen-only; TPS/non-TPS hits)."},
    "maxid_gen_vs_train":   {"default": True,  "branch": "sequence",  "description": "Max sequence identity of each gen seq vs the train set."},
    "mindist_gen_vs_train": {"default": True,  "branch": "sequence",  "description": "Min ESM-embedding distance of gen vs train (needs esm)."},
    "esmfold":              {"default": False, "branch": "producer",  "description": "ESMFold structure PRODUCER (both clusters): folds the gen FASTA into a structs dir + PAE. Opt-in via --fold esmfold (auto-enabled then); not a default-on metric."},
    "alphafold3":           {"default": False, "branch": "producer",  "description": "AlphaFold3 structure PRODUCER (Aurum-only): per-sequence fan-out (one AF3 job each) into <gen>_af3/structs + af_output, then a PAE-extraction step. Opt-in via --fold alphafold3 (auto-enabled then); not a default-on metric."},
    "plddt":                {"default": True,  "branch": "structure", "description": "AlphaFold/ESMFold pLDDT folding confidence."},
    "motif_struct":         {"default": True,  "branch": "structure", "description": "Structural distance between the two metal-binding motifs."},
    "active_site_geom":     {"default": True,  "branch": "structure", "description": "Active-site carboxylate-cage geometry (apo-robust)."},
    "aromatic_lining":      {"default": True,  "branch": "structure", "description": "Aromatic / cation-pi pocket lining (carbocation-stabilization proxy)."},
    "diphosphate_sensor":   {"default": True,  "branch": "structure", "description": "Diphosphate-sensor basic residues (Arg/Lys + RY pair) at the metal site."},
    "ion_site":             {"default": True,  "branch": "structure", "description": "Ion-placement check: do AF3 co-folded Mg/Mn ions land in the carboxylate cage? Only carries signal for AF3 holo folds (--af3_cofold mg*); apo/ESMFold report n_ions_modelled=0. Gated on a holo co-fold / --no_holo_tools."},
    "substrate_positioning":{"default": True,  "branch": "structure", "description": "Substrate positioning: is the AF3 co-folded prenyl-PP substrate poised in the catalytic cage (diphosphate->FARM/Mg, reactive C1->cage)? Auto-detects the ligand per design (--af3_cofold mg_<sub>|mg_ee); apo / no substrate -> NaN. Gated on a holo co-fold / --no_holo_tools."},
    "radius_of_gyration":   {"default": True,  "branch": "structure", "description": "Radius of gyration / compactness over Ca atoms."},
    "pocket_descriptors":   {"default": True,  "branch": "structure", "description": "Active-site pocket descriptors (fpocket volume/hydrophobicity/enclosure + P2Rank ligandability cross-check)."},
    "domain_composition":   {"default": True,  "branch": "structure", "description": "TPS structural-domain composition (EE CPU detector)."},
    "aggregation":          {"default": True,  "branch": "structure", "description": "Aggrescan3D structure-based aggregation propensity."},
    "foldseek_swissprot":   {"default": True,  "branch": "structure", "description": "Foldseek search vs AlphaFold-Swiss-Prot (TPS/non-TPS hits)."},
    "structural_identity":  {"default": True,  "branch": "structure", "description": "Foldseek structural identity to nearest known TPS (needs --known_structs_dir)."},
    "proteinmpnn":          {"default": True,  "branch": "structure", "description": "ProteinMPNN sequence-likelihood (NLL) of the design's own sequence given its fold."},
    "global_confidence":    {"default": True,  "branch": "structure", "description": "Global fold confidence (pTM/iPTM) from the saved PAE npz (needs --pae_dir)."},
    "interdomain_pae":      {"default": True,  "branch": "structure", "description": "Mean/max inter-domain PAE between TPS domains (needs --pae_dir; EE domain ranges)."},
    "self_consistency":     {"default": False, "branch": "structure", "description": "HEAVY scRMSD self-consistency (ProteinMPNN -> ESMFold refold -> RMSD). Opt-in."},
    "knn":                  {"default": True,  "branch": "structure", "description": "k-NN coarse-label transfer: ensembled vote over the sequence/embedding/structural top-k neighbours -> predicted coarse label + confidence (needs --train_path + --structs_dir + --known_structs_dir)."},
    "sdr_divergence":       {"default": True,  "branch": "structure", "description": "SDR specificity-divergence: flags designs globally close to a known TPS but divergent at the specificity-determining active-site residues (needs --structs_dir + --known_structs_dir)."},
    "domain_structural_identity": {"default": True, "branch": "structure", "description": "Domain-level structural identity: EE detects each design's TPS domains, then foldseek-aligns them to the known martsDB reference domains (per-domain-type best TM-score/lddt; n_detected_domains)."},
    "substrate_class":      {"default": True,  "branch": "structure", "description": "Substrate-class combiner: fuses the SUBSTRATE k-NN vote (3 spaces) with the pocket-volume size band + EnzymeExplorer per-substrate signal -> predicted substrate (GPP/FPP/GGPP/...) + agreement (needs --train_path + --structs_dir + --known_structs_dir)."},
    "plots":                {"default": True,  "branch": "sequence",  "description": "Aggregator: merges all enabled metrics into plots. Effectively always on unless excluded."},
}


def load_tools_catalog() -> Dict[str, dict]:
    """Load the tools catalog from JSON, falling back to DEFAULT_TOOLS.

    JSON keys that exist in DEFAULT_TOOLS override the built-in entry; any
    DEFAULT_TOOLS key missing from the JSON is filled from the fallback (so a stale
    config that simply lacks a freshly-added tool still runs that tool). Keys
    starting with '_' (e.g. '_comment') are ignored.
    """
    catalog = {k: dict(v) for k, v in DEFAULT_TOOLS.items()}
    try:
        with open(TOOLS_CONFIG) as fh:
            data = json.load(fh)
    except (OSError, ValueError) as e:
        print(f"[note] tools config not used ({TOOLS_CONFIG}: {e}); using built-in defaults.")
        return catalog
    for key, spec in data.items():
        if key.startswith("_"):
            continue
        if not isinstance(spec, dict):
            continue
        entry = catalog.get(key, {"default": True, "branch": "sequence", "description": ""})
        entry.update({k: spec[k] for k in ("default", "branch", "description") if k in spec})
        catalog[key] = entry
    return catalog


def knn_calibrated_spaces(calibration_path: str) -> set:
    """Return the similarity spaces the k-NN calibration JSON actually covers.

    `transfer_labels` indexes calibration["spaces"][space] for EVERY space whose
    top-k CSV is passed, so passing an UNCALIBRATED space (e.g. 'sequence' while only
    embedding+structural are calibrated) raises KeyError. The orchestrator therefore
    only forwards the topk flags for the calibrated spaces. Returns an empty set if the
    file is unreadable (the caller then [note]-skips knn).
    """
    try:
        with open(calibration_path) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return set()
    spaces = data.get("spaces")
    return set(spaces) if isinstance(spaces, dict) else set()


def load_top_k() -> int:
    """Resolve the pipeline-wide top-k from JSON "_settings.top_k", else DEFAULT_TOP_K.

    The tools-catalog loader ignores '_'-prefixed keys, so the top-k setting is read
    here explicitly from "_settings": {"top_k": N}. Missing/unreadable/malformed ->
    the built-in DEFAULT_TOP_K fallback.
    """
    try:
        with open(TOOLS_CONFIG) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return DEFAULT_TOP_K
    settings = data.get("_settings")
    if isinstance(settings, dict):
        val = settings.get("top_k")
        if isinstance(val, int) and val >= 1:
            return val
    return DEFAULT_TOP_K


def _parse_keys(csv: Optional[str]) -> List[str]:
    return [t.strip() for t in csv.split(",") if t.strip()] if csv else []


def resolve_enabled_tools(catalog: Dict[str, dict], args) -> set:
    """Compute the set of ENABLED tool keys from config defaults + CLI overrides.

    Precedence: config defaults -> --only (restrict) -> --include (add) -> --exclude
    (remove). `plots` is auto-kept on unless it is explicitly excluded (or --only is
    given without it AND it is not later excluded -> we re-add it).
    Unknown tool keys raise SystemExit listing the valid keys.
    """
    valid = set(catalog)
    only = _parse_keys(args.only)
    include = _parse_keys(args.include)
    exclude = _parse_keys(args.exclude)
    if args.self_consistency:                       # back-compat alias
        include.append("self_consistency")

    unknown = [k for k in (only + include + exclude) if k not in valid]
    if unknown:
        raise SystemExit(
            f"[FAIL] unknown tool key(s): {', '.join(sorted(set(unknown)))}\n"
            f"       valid keys: {', '.join(sorted(valid))}")

    if only:
        enabled = set(only)
        enabled.add("plots")                        # aggregator stays on unless excluded
    else:
        enabled = {k for k, v in catalog.items() if v.get("default", False)}
    enabled |= set(include)
    # --fold turns on the (opt-in) producer tool so its Step(s) survive the tool filter.
    # Ignored when --structs_dir is given (structures already exist; no folding).
    fold = getattr(args, "fold", None)
    if fold and not args.structs_dir:
        enabled.add("esmfold" if fold == "esmfold" else "alphafold3")
    enabled -= set(exclude)
    return enabled


# --------------------------------------------------------------------------- #
# Step model + engine                                                         #
# --------------------------------------------------------------------------- #
@dataclass
class Step:
    name: str
    job: str                       # job-script basename under scripts/<cluster>/jobs/
    args: List[str]
    output: str                    # file OR dir whose existence means "already done"
    tool: str = ""                 # tool KEY this step belongs to (config-driven on/off)
    deps: List[str] = field(default_factory=list)        # HARD deps: skip step if unmet
    soft_deps: List[str] = field(default_factory=list)   # wait-on-if-submitted; never gates
    requires_cap: Optional[str] = None
    driver: bool = False           # login-node FAN-OUT driver: run args directly (no
                                   # sbatch), capturing the N submitted job ids it prints.
                                   # Downstream deps on this step wait on ALL N ids.


class Engine:
    def __init__(self, cluster: str, log_dir: str, dry_run: bool = False):
        self.cfg = CLUSTERS[cluster]
        self.cluster = cluster
        self.log_dir = log_dir
        self.dry_run = dry_run
        self.jobs_dir = os.path.join(REPO, "scripts", cluster, "jobs")
        # step -> list of SLURM ids (normal steps -> [id]; fan-out driver -> [id1, id2, …];
        # skipped/output-exists steps -> absent, so downstream deps add nothing and don't wait).
        self.job_ids: Dict[str, List[str]] = {}
        self.satisfied: set = set()         # steps whose output exists or were submitted

    def _exists(self, path: str) -> bool:
        if os.path.isdir(path):
            return any(os.scandir(path))    # dir outputs count only if non-empty
        return os.path.isfile(path)

    def run(self, steps: List[Step]) -> None:
        for s in steps:
            if s.requires_cap and s.requires_cap not in self.cfg["caps"]:
                print(f"[cap ] {s.name}: skipped ('{s.requires_cap}' unavailable on {self.cluster})")
                continue
            missing = [d for d in s.deps if d not in self.satisfied]
            if missing:
                print(f"[dep ] {s.name}: skipped (unsatisfied deps: {', '.join(missing)})")
                continue
            if self._exists(s.output):
                print(f"[skip] {s.name}: output exists ({os.path.relpath(s.output, REPO)})")
                self.satisfied.add(s.name)
                continue

            # Dep ids: a dep that fanned out contributes ALL its job ids (flatten).
            dep_ids = [jid for d in (s.deps + s.soft_deps)
                       for jid in self.job_ids.get(d, [])]

            # Fan-out driver: NOT an sbatch job — run the command directly on the login
            # node (the orchestrator already runs there). It internally submits N jobs and
            # prints their ids; capture all so downstream steps afterok-wait on every one.
            if s.driver:
                if self.dry_run:
                    wait = f" [after {':'.join(dep_ids)}]" if dep_ids else ""
                    print(f"[dry ] {s.name} (fan-out driver){wait}: {' '.join(s.args)}")
                    self.job_ids[s.name] = [f"<{s.name}:N>"]
                    self.satisfied.add(s.name)
                    continue
                proc = subprocess.run(s.args, capture_output=True, text=True)
                out = (proc.stdout or "") + (proc.stderr or "")
                if proc.returncode != 0:
                    raise SystemExit(f"[FAIL] {s.name}: fan-out driver failed\n{out}")
                jids = self.cfg["fanout_ids"](out)
                print(f"[fan ] {s.name}: submitted {len(jids)} fold job(s)"
                      f"{' ' + ':'.join(jids) if jids else ' (none — all structures present)'}")
                self.job_ids[s.name] = jids
                self.satisfied.add(s.name)
                continue

            job_path = os.path.join(self.jobs_dir, s.job)
            if not os.path.isfile(job_path):
                print(f"[MISS] {s.name}: job script not found ({job_path}) -- skipping")
                continue

            cmd = ([self.cfg["submit"]] + self.cfg["dep"](dep_ids)
                   + self.cfg["log"](self.log_dir) + [job_path] + s.args)
            if self.dry_run:
                wait = f" [after {':'.join(dep_ids)}]" if dep_ids else ""
                print(f"[dry ] {s.name}{wait}: {' '.join(cmd)}")
                self.job_ids[s.name] = [f"<{s.name}>"]
                self.satisfied.add(s.name)
                continue
            proc = subprocess.run(cmd, capture_output=True, text=True)
            out = (proc.stdout or "") + (proc.stderr or "")
            try:
                jid = self.cfg["jobid"](out)
            except RuntimeError as e:
                raise SystemExit(f"[FAIL] {s.name}: submission failed\n{out}") from e
            dep_note = f" (after {':'.join(dep_ids)})" if dep_ids else ""
            print(f"[sub ] {s.name}: job {jid}{dep_note}")
            self.job_ids[s.name] = [jid]
            self.satisfied.add(s.name)


# --------------------------------------------------------------------------- #
# Pipeline definition                                                         #
# --------------------------------------------------------------------------- #
def build_steps(args, enabled: set) -> List[Step]:
    """Build the full candidate Step list, then keep only steps whose tool key is
    enabled. `plots` (the aggregator) soft-depends only on the metric steps that
    survived filtering, so it never waits on a disabled/absent tool. Hard-dep
    semantics are unchanged: if a step's enabled but its hard dep (e.g. esm for
    mindist_*) was filtered out, the Engine emits the existing '[dep ]' skip
    message rather than auto-pulling the dep -- see resolve_enabled_tools docs."""
    gen = args.fasta_path
    train = args.train_path
    steps: List[Step] = []

    structs = args.structs_dir
    known_structs = args.known_structs_dir
    pae_dir = args.pae_dir

    # Fold producer: when --fold is given without a pre-supplied --structs_dir, the
    # producer makes the structures (+ PAE) the structure branch then consumes.
    #   esmfold     -> one whole-FASTA job; dirs mirror run_esmfold.sh's defaults.
    #   alphafold3  -> a per-sequence FAN-OUT (one AF3 job each, Aurum-only) under a
    #                  <gen>_af3/ work dir (structs/ + af_output/), then an extract_pae step.
    fold_mode = getattr(args, "fold", None)
    do_fold = bool(fold_mode) and not args.structs_dir
    af3_work = None
    if do_fold and fold_mode == "esmfold":
        structs = _base(gen) + "_esmfold_structs"
        if pae_dir is None and not args.no_fold_pae:
            pae_dir = structs + "_pae"
    elif do_fold and fold_mode == "alphafold3":
        af3_work = _base(gen) + "_af3"
        structs = os.path.join(af3_work, "structs")  # where alphafold.sh writes <ID>.pdb
        if pae_dir is None and not args.no_fold_pae:
            pae_dir = os.path.join(af3_work, "pae")

    if args.af3_cofold == "mg_ee" and not args.enzymeexplorer_csv:
        raise SystemExit(
            "[FAIL] --af3_cofold mg_ee needs --enzymeexplorer_csv (the EnzymeExplorer seq-only CSV): the "
            "AF3 fan-out runs on the login node and cannot wait on the in-pipeline ee_seq job. "
            "Run enzyme_explorer_sequence_only first, then pass its CSV via --enzymeexplorer_csv.")
    # Holo (co-fold-dependent) tools run iff we co-folded a holo active site (--af3_cofold !=
    # none) OR structures were supplied externally (which may be holo) -- and never when
    # --no_holo_tools is set. So --af3_cofold none (default) cleanly turns co-folding AND the
    # downstream holo tools (ion_site_check, substrate_positioning) off.
    run_holo = (not args.no_holo_tools) and (args.af3_cofold != "none" or bool(args.structs_dir))

    if args.train_structs_dir:
        print("[warn] --train_structs_dir is not used by the orchestrator yet "
              "(EnzymeExplorer-with-structures is not ported).")

    datasets = [("gen", gen)] + ([("train", train)] if train else [])

    # The k-NN + SDR consumers read the gen-vs-train / gen-vs-known *_topk.csv CSVs
    # emitted by the three feeders: local_sequence_search (sequence, MMseqs2 — fast),
    # min_embedding_distance (embedding), structural_identity (structural). Only ask the
    # feeders for top-k when a consumer is actually enabled (no wasted top-k otherwise).
    # topk_args is appended to the three gen-vs-{train,known} feeder Steps below.
    # max_sequence_identity is NO LONGER a feeder — it stays as the plain GLOBAL novelty
    # metric (no --top_k); the fast local search supplies the sequence-space top-k.
    topk_consumers = {"knn", "sdr_divergence", "substrate_class"}
    want_topk = bool(topk_consumers & enabled)
    topk_args = ["--top_k", str(args.top_k)] if want_topk else []

    for tag, fa in datasets:
        is_train = tag == "train"
        steps.append(Step(f"motif_{tag}", "motif_search.sh", ["--fasta_path", fa],
                          out_motif(fa), tool="motif"))
        steps.append(Step(f"motif_pair_{tag}", "motif_pair_distance.sh",
                          ["--fasta_path", fa], out_motif_pair(fa), tool="motif_pair"))
        steps.append(Step(f"esm_{tag}", "esm_embedding.sh", ["--fasta_path", fa],
                          out_esm(fa), tool="esm"))
        steps.append(Step(f"esm_ppl_{tag}", "esm_pseudo_perplexity.sh",
                          ["--fasta_path", fa], out_esm_ppl(fa), tool="esm_ppl"))
        steps.append(Step(f"maxid_self_{tag}", "max_sequence_identity.sh",
                          ["--fasta_path", fa] + (["--train"] if is_train else []),
                          out_maxid_self(fa), tool="maxid_self"))
        # Fast LOCAL (MMseqs2) sequence search within the dataset (self mode -> no
        # --train_path; the tool excludes each query's self-hit). Within-set novelty
        # counterpart of maxid_self. Explicit --save_path -> _self.csv so it does not
        # collide with the gen-vs-train step's default <fa>_local_sequence_search.csv.
        steps.append(Step(f"local_search_{tag}", "local_sequence_search.sh",
                          ["--fasta_path", fa, "--save_path", out_local_search_self(fa)],
                          out_local_search_self(fa), tool="local_sequence_search"))
        steps.append(Step(f"mindist_self_{tag}", "min_embedding_distance.sh",
                          ["--embeddings_path", out_esm(fa)] + (["--train"] if is_train else []),
                          out_mindist_self(fa), tool="mindist_self", deps=[f"esm_{tag}"]))
        steps.append(Step(f"soluprot_{tag}", "soluprot.sh", ["--fasta_path", fa],
                          out_soluprot(fa), tool="soluprot"))
        steps.append(Step(f"ee_seq_{tag}", "enzyme_explorer_sequence_only.sh",
                          ["--fasta_path", fa], out_ee_seq(fa), tool="ee_seq"))

    # Broad "what-else-is-it-like" sequence search vs Swiss-Prot (annotated -> each hit
    # labeled TPS/non-TPS). Gen-only: it evaluates the generated designs (real train TPS
    # are trivially TPS hits).
    steps.append(Step("swissprot_search_gen", "swissprot_search.sh",
                      ["--fasta_path", gen], out_swissprot_search(gen), tool="swissprot_search"))

    if train:
        train_emb = args.train_embeddings_path or out_esm(train)
        # max_sequence_identity stays as the plain GLOBAL novelty metric (no --top_k):
        # the fast local_sequence_search now feeds the sequence-space top-k instead.
        steps.append(Step("maxid_gen_vs_train", "max_sequence_identity.sh",
                          ["--fasta_path", gen, "--train_path", train], out_maxid(gen),
                          tool="maxid_gen_vs_train"))
        # Fast LOCAL (MMseqs2) gen-vs-train search. Its _topk.csv (query_id,rank,
        # neighbour_id,score=identity%) is the SEQUENCE-space feeder for knn + sdr.
        steps.append(Step("local_search_gen_vs_train", "local_sequence_search.sh",
                          ["--fasta_path", gen, "--train_path", train] + topk_args,
                          out_local_search(gen), tool="local_sequence_search"))
        steps.append(Step("mindist_gen_vs_train", "min_embedding_distance.sh",
                          ["--embeddings_path", out_esm(gen), "--train_embeddings_path", train_emb]
                          + topk_args,
                          out_mindist(gen), tool="mindist_gen_vs_train", deps=["esm_gen", "esm_train"]))

    # Structure branch: pLDDT (folding confidence) + foldseek structural identity to the
    # nearest known TPS. Both are keyed by the structures dir, run if --structs_dir is
    # given OR --fold produced one. With --fold the structures don't exist out-of-band, so
    # every structure step waits on the producer (attached by the fixup at the end of the
    # block); without --fold they have no SLURM dep and plots picks them up via the snapshot.
    if structs:
        struct_block_start = len(steps)
        fold_producer = None          # step name the structure branch waits on
        pae_producer = None           # step name the PAE-consumers wait on (AF3: a later step)
        if do_fold and fold_mode == "esmfold":
            # ESMFold producer: folds the gen FASTA into `structs` (+ PAE into pae_dir).
            # output=structs (a dir) -> the Engine skips folding if it already exists
            # and is non-empty (idempotent re-runs). PAE is written by the SAME job.
            esmfold_args = ["--fasta_path", gen, "--save_dir", structs]
            esmfold_args += (["--pae_dir", pae_dir] if pae_dir else ["--no-save_pae"])
            steps.append(Step("esmfold_gen", "esmfold.sh", esmfold_args, structs,
                              tool="esmfold"))
            fold_producer = pae_producer = "esmfold_gen"
        elif do_fold and fold_mode == "alphafold3":
            # AF3 FAN-OUT producer (Aurum-only). A login-node driver runs the existing
            # run_alphafold_jobs.py (one AF3 job per sequence -> structs/ + af_output/) and
            # prints the N job ids; the Engine captures them so the structure branch waits
            # on all N. A sentinel output (never created) -> the driver always runs and
            # run_alphafold_jobs' own --skip_existing handles per-design idempotency.
            fanout_cmd = ["bash", os.path.join(REPO, "scripts", "run_alphafold_fanout.sh"),
                          "--cluster", args.cluster, "--fasta_path", gen,
                          "--working_directory", af3_work,
                          "--cofold", args.af3_cofold]
            if args.af3_cofold == "mg_ee":
                # Per-design EE substrate: the login-node driver needs the EE CSV up front
                # (it cannot afterok-wait on the in-pipeline ee_seq SLURM job).
                fanout_cmd += ["--enzymeexplorer_csv", args.enzymeexplorer_csv]
            fanout_cmd += ["--model_seeds"] + [str(s) for s in args.af3_model_seeds]
            steps.append(Step("af3_fold_gen", "", fanout_cmd,
                              os.path.join(af3_work, ".__never__"),
                              tool="alphafold3", driver=True, requires_cap="alphafold"))
            fold_producer = "af3_fold_gen"
            if pae_dir:
                # Extract <ID>_pae.npz from the af_output tree once ALL fold jobs finish.
                steps.append(Step("af3_pae_gen", "extract_pae.sh",
                                  ["--structs_dir", structs, "--pae_dir", pae_dir],
                                  pae_dir, tool="alphafold3", deps=["af3_fold_gen"],
                                  requires_cap="alphafold"))
                pae_producer = "af3_pae_gen"
        steps.append(Step("plddt_gen", "plddt.sh", ["--structs_dir", structs],
                          out_plddt(structs), tool="plddt"))
        steps.append(Step("motif_struct_gen", "motif_structural_distance.sh",
                          ["--structs_dir", structs], out_motif_struct(structs), tool="motif_struct"))
        # Active-site geometry: carboxylate-cage convergence (apo-robust). The
        # constellation-RMSD columns need reference templates present in the structs
        # dir (pass --templates <ids>); omitted here until a curated reference set is wired.
        steps.append(Step("active_site_geom_gen", "active_site_geometry.sh",
                          ["--structs_dir", structs], out_active_site_geom(structs),
                          tool="active_site_geom"))
        # Aromatic / cation-π pocket lining (carbocation-stabilization proxy).
        steps.append(Step("aromatic_lining_gen", "aromatic_lining.sh",
                          ["--structs_dir", structs], out_aromatic_lining(structs),
                          tool="aromatic_lining"))
        # Diphosphate-sensor basic residues (Arg/Lys + RY pair) at the metal site.
        steps.append(Step("diphosphate_sensor_gen", "diphosphate_sensor.sh",
                          ["--structs_dir", structs], out_diphosphate_sensor(structs),
                          tool="diphosphate_sensor"))
        # Holo (co-fold-dependent) tools -- only meaningful when ions/substrate are modelled.
        # Gated on run_holo (--af3_cofold != none, or external structs; off via --no_holo_tools)
        # so that apo runs don't queue no-op jobs.
        if run_holo:
            # Ion-placement check: do the AF3 co-folded Mg/Mn ions (--af3_cofold mg*) actually
            # land in the carboxylate cage? apo/ESMFold structures report n_ions_modelled=0.
            steps.append(Step("ion_site_check_gen", "ion_site_check.sh",
                              ["--structs_dir", structs], out_ion_site_check(structs),
                              tool="ion_site"))
            # Substrate positioning: is the co-folded prenyl-PP substrate poised in the cage?
            # Auto-detects the ligand per design (works for forced mg_<sub> and per-design
            # mg_ee alike); reports NaN when no substrate ligand is present.
            steps.append(Step("substrate_positioning_gen", "substrate_positioning.sh",
                              ["--structs_dir", structs], out_substrate_positioning(structs),
                              tool="substrate_positioning"))
        # Radius of gyration / compactness: raw geometric shape numbers (Rg, asphericity,
        # principal radii) over the Cα atoms; no expected-Rg band (compared downstream).
        steps.append(Step("radius_of_gyration_gen", "radius_of_gyration.sh",
                          ["--structs_dir", structs], out_radius_of_gyration(structs),
                          tool="radius_of_gyration"))
        # Active-site pocket descriptors: fpocket geometry (catalytic-pocket volume,
        # hydrophobicity, enclosure) + P2Rank ML ligandability, anchored on the
        # carboxylate-cage metal point. Raw numbers; band from the reference-stats pipeline.
        steps.append(Step("pocket_descriptors_gen", "pocket_descriptors.sh",
                          ["--structs_dir", structs], out_pocket_descriptors(structs),
                          tool="pocket_descriptors"))
        # TPS structural-domain composition via EnzymeExplorer's CPU domain detector.
        steps.append(Step("domain_composition_gen", "domain_composition.sh",
                          ["--structs_dir", structs], out_domain_composition(structs),
                          tool="domain_composition"))
        # Domain-level structural identity: detect each design's TPS domains (EE) and
        # foldseek-align them to the known martsDB reference domains. Reference domain
        # root defaults inside run_domain_structural_identity.sh to EE's curated set, so
        # the Step needs only --structs_dir (gen-only; no --known_structs_dir).
        steps.append(Step("domain_structural_identity_gen", "domain_structural_identity.sh",
                          ["--structs_dir", structs], out_domain_structural_identity(structs),
                          tool="domain_structural_identity"))
        # Aggrescan3D structure-based aggregation propensity (expressibility signal).
        steps.append(Step("aggregation_gen", "aggregation.sh",
                          ["--structs_dir", structs], out_aggregation(structs), tool="aggregation"))
        # Broad structural search vs AlphaFold-Swiss-Prot (annotated -> hit TPS/non-TPS).
        steps.append(Step("foldseek_swissprot_gen", "foldseek_swissprot_search.sh",
                          ["--structs_dir", structs], out_foldseek_swissprot(structs),
                          tool="foldseek_swissprot"))
        # ProteinMPNN sequence-likelihood (NLL) of the design's own sequence given its fold.
        steps.append(Step("proteinmpnn_gen", "proteinmpnn_score.sh",
                          ["--structs_dir", structs], out_proteinmpnn(structs), tool="proteinmpnn"))
        # Self-consistency scRMSD (ProteinMPNN -> ESMFold refold -> RMSD). HEAVY
        # (~1-2.5 min/structure x num_seqs GPU), so opt-in (config default off;
        # enable via --include self_consistency or the --self_consistency alias).
        steps.append(Step("self_consistency_gen", "self_consistency.sh",
                          ["--structs_dir", structs, "--num_seqs", str(args.self_consistency_num_seqs)],
                          out_self_consistency(structs), tool="self_consistency"))
        if known_structs:
            steps.append(Step("structural_identity_gen", "structural_identity.sh",
                              ["--structs_dir", structs, "--known_structs_dir", known_structs]
                              + topk_args,
                              out_structural_identity(structs), tool="structural_identity"))
        else:
            print("[note] --structs_dir without --known_structs_dir: skipping "
                  "structural_identity (needs a known-TPS reference structure dir).")
        # PAE-consuming metrics need the per-structure PAE matrices saved at fold time
        # (<ID>_pae.npz, produced by ESMFold / the AF3 extract_pae step). Gated on
        # --pae_dir, like structural_identity is gated on --known_structs_dir.
        if pae_dir:
            steps.append(Step("global_confidence_gen", "global_confidence.sh",
                              ["--structs_dir", structs, "--pae_dir", pae_dir],
                              out_global_confidence(structs), tool="global_confidence"))
            steps.append(Step("interdomain_pae_gen", "interdomain_pae.sh",
                              ["--structs_dir", structs, "--pae_dir", pae_dir],
                              out_interdomain_pae(structs), tool="interdomain_pae"))
        else:
            print("[note] --structs_dir without --pae_dir: skipping global_confidence + "
                  "interdomain_pae (need <ID>_pae.npz from a PAE-saving fold).")

        # With --fold the structures are produced in-pipeline, so every structure step in
        # this block must wait on the producer (the producer steps themselves excepted).
        # PAE-consumers wait on the PAE producer specifically (for AF3 that's the separate
        # af3_pae_gen extraction step, which itself waits on all fold jobs; for ESMFold the
        # fold job writes the PAE so it's the same step). The downstream knn/sdr/substrate
        # consumers need nothing added — they hard-dep on structural_identity_gen (which
        # gets the producer dep here) and SLURM afterok is transitive in effect.
        if do_fold:
            pae_consumers = {"global_confidence_gen", "interdomain_pae_gen"}
            producers = {p for p in (fold_producer, pae_producer) if p}
            for s in steps[struct_block_start:]:
                if s.name in producers:
                    continue
                dep = pae_producer if (s.name in pae_consumers and pae_producer) else fold_producer
                if dep and dep not in s.deps:
                    s.deps = [dep] + list(s.deps)

    # --- Dependency-heavy consumers: k-NN label transfer + SDR divergence --------- #
    # Both read the three feeders' *_topk.csv. They are gated on the feeders actually
    # existing (train -> sequence/embedding feeders; structs+known_structs -> structural
    # feeder) and skip cleanly with a [note] otherwise, like structural_identity does
    # without --known_structs_dir.
    #
    # k-NN: ensembled coarse-label transfer over the three similarity spaces. predict
    # mode needs ALL THREE feeders present, plus the label_file + calibration artifacts.
    # The committed calibration covers all three spaces (sequence via the fast MMseqs2
    # local_sequence_search). knn_calibrated_spaces() forwards only the calibrated spaces,
    # so a 2-space calibration would still work (sequence space abstained).
    if "knn" in enabled:
        knn_ok = bool(train and structs and known_structs)
        if not knn_ok:
            print("[note] knn: skipping label transfer (needs --train_path AND "
                  "--structs_dir AND --known_structs_dir for the three feeders).")
        elif not os.path.isfile(args.knn_label_file):
            print(f"[note] knn: skipping (label file not found: {args.knn_label_file}).")
        elif not os.path.isfile(args.knn_calibration):
            print(f"[note] knn: skipping (calibration not found: {args.knn_calibration}).")
        else:
            # Only forward the topk for spaces the calibration actually covers:
            # transfer_labels indexes calibration["spaces"][space] for every passed
            # space, so an uncalibrated one (e.g. sequence while only embedding+
            # structural are calibrated) would KeyError. space -> (flag, topk, feeder).
            cal_spaces = knn_calibrated_spaces(args.knn_calibration)
            space_specs = {
                "sequence":   ("--sequence_topk",   out_local_search_topk(gen),
                               "local_search_gen_vs_train"),
                "embedding":  ("--embedding_topk",  out_mindist_topk(out_esm(gen)),
                               "mindist_gen_vs_train"),
                "structural": ("--structural_topk", out_structural_identity_topk(structs),
                               "structural_identity_gen"),
            }
            knn_args = ["predict"]
            knn_deps: List[str] = []
            for space in ("sequence", "embedding", "structural"):
                if space not in cal_spaces:
                    continue
                flag, topk_path, feeder = space_specs[space]
                knn_args += [flag, topk_path]
                knn_deps.append(feeder)
            if not knn_deps:
                print(f"[note] knn: skipping (calibration covers no known space: "
                      f"{sorted(cal_spaces)}).")
            else:
                if "sequence" not in cal_spaces:
                    print("[note] knn: sequence space not yet calibrated -> using "
                          f"{sorted(s for s in cal_spaces if s in space_specs)} only.")
                knn_args += ["--label_file", args.knn_label_file,
                             "--calibration", args.knn_calibration,
                             "--output", out_knn(gen)]
                steps.append(Step("knn_gen", "knn_label_transfer.sh", knn_args,
                                  out_knn(gen), tool="knn", deps=knn_deps))

    # SDR divergence: needs the structural feeder (structs + known_structs); the sequence
    # feeder (--sequence_topk) is an optional fallback added only when train is present.
    if "sdr_divergence" in enabled:
        if not (structs and known_structs):
            print("[note] sdr_divergence: skipping (needs --structs_dir AND "
                  "--known_structs_dir).")
        else:
            sdr_args = ["--structs_dir", structs, "--known_structs_dir", known_structs,
                        "--structural_topk", out_structural_identity_topk(structs)]
            sdr_deps = ["structural_identity_gen"]
            if train:
                sdr_args += ["--sequence_topk", out_local_search_topk(gen)]
                sdr_deps.append("local_search_gen_vs_train")
            sdr_args += ["--save_path", out_sdr_divergence(structs)]
            steps.append(Step("sdr_divergence_gen", "sdr_divergence.sh", sdr_args,
                              out_sdr_divergence(structs), tool="sdr_divergence", deps=sdr_deps))

    # Substrate-class combiner: same three feeders as k-NN (run with the SUBSTRATE label
    # file + calibration), cross-checked against pocket_descriptors volume + the EE
    # per-substrate sequence signal. Gated like k-NN (needs all three spaces -> train +
    # structs + known_structs). Forwards only the spaces the substrate calibration covers.
    if "substrate_class" in enabled:
        sc_ok = bool(train and structs and known_structs)
        if not sc_ok:
            print("[note] substrate_class: skipping (needs --train_path AND --structs_dir "
                  "AND --known_structs_dir for the three feeders).")
        elif not os.path.isfile(args.substrate_label_file):
            print(f"[note] substrate_class: skipping (label file not found: "
                  f"{args.substrate_label_file}).")
        elif not os.path.isfile(args.substrate_calibration):
            print(f"[note] substrate_class: skipping (calibration not found: "
                  f"{args.substrate_calibration}).")
        else:
            cal_spaces = knn_calibrated_spaces(args.substrate_calibration)
            sc_space_specs = {
                "sequence":   ("--sequence_topk",   out_local_search_topk(gen),
                               "local_search_gen_vs_train"),
                "embedding":  ("--embedding_topk",  out_mindist_topk(out_esm(gen)),
                               "mindist_gen_vs_train"),
                "structural": ("--structural_topk", out_structural_identity_topk(structs),
                               "structural_identity_gen"),
            }
            sc_args: List[str] = []
            sc_deps: List[str] = []
            for space in ("sequence", "embedding", "structural"):
                if space not in cal_spaces:
                    continue
                flag, topk_path, feeder = sc_space_specs[space]
                sc_args += [flag, topk_path]
                sc_deps.append(feeder)
            if not sc_deps:
                print(f"[note] substrate_class: skipping (calibration covers no known "
                      f"space: {sorted(cal_spaces)}).")
            else:
                # pocket-volume cross-check (structs) + EE per-substrate signal (gen fasta).
                sc_args += ["--pocket_csv", out_pocket_descriptors(structs),
                            "--enzymeexplorer_csv", out_ee_seq(gen)]
                sc_deps += ["pocket_descriptors_gen", "ee_seq_gen"]
                sc_args += ["--top_k", str(args.top_k),
                            "--label_file", args.substrate_label_file,
                            "--calibration", args.substrate_calibration,
                            "--output", out_substrate_class(gen)]
                steps.append(Step("substrate_class_gen", "substrate_class.sh", sc_args,
                                  out_substrate_class(gen), tool="substrate_class", deps=sc_deps))

    # Filter to enabled tools BEFORE wiring plots, so plots only soft-depends on the
    # metric steps that actually survived (it never waits on a disabled tool).
    steps = [s for s in steps if s.tool in enabled]

    # plots: depend SOFTLY on every (enabled) metric (wait on whatever was submitted;
    # still run if some metrics were skipped/absent — the plots tool skips
    # missing-input targets). Only added if 'plots' itself is enabled.
    if "plots" in enabled:
        metric_steps = [s.name for s in steps]
        plot_dir = os.path.join(os.path.dirname(os.path.abspath(gen)), "plots")
        if train:
            plot_args = ["--fasta_paths", train, gen, "--data_names", "train", "generated",
                         "--data_colors", args.data_colors[0], args.data_colors[1],
                         "--save_dir", plot_dir]
        else:
            plot_args = ["--fasta_paths", gen, "--data_names", "generated",
                         "--data_colors", args.data_colors[1], "--save_dir", plot_dir]
        # Sentinel output that never exists -> plots always (re)runs after the metrics.
        steps.append(Step("plots", "plots.sh", plot_args, os.path.join(plot_dir, ".__never__"),
                          tool="plots", soft_deps=metric_steps))
    return steps


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cluster", required=True, choices=sorted(CLUSTERS))
    p.add_argument("--account", default=None,
                   help="SLURM account for sbatch (exported as $SBATCH_ACCOUNT). Default: an "
                        "existing $SBATCH_ACCOUNT, else scripts/<cluster>/config.sh (which sources "
                        "the per-install scripts/<cluster>/config.local.sh).")
    p.add_argument("--fasta_path", required=True, help="Generated sequences FASTA.")
    p.add_argument("--train_path", default=None, help="Reference/train FASTA (enables comparisons).")
    p.add_argument("--structs_dir", default=None,
                   help="Generated structures dir (AF3 af_output or flat .pdb/.cif) -> "
                        "enables pLDDT and (with --known_structs_dir) structural-identity steps.")
    p.add_argument("--fold", default=None, choices=["esmfold", "alphafold3"],
                   help="Fold the generated sequences FIRST, then run the structure branch "
                        "on the produced structures (no need to pre-supply --structs_dir). "
                        "'esmfold' runs ESMFold (both clusters) into <gen>_esmfold_structs/ "
                        "+ a sibling _pae/. 'alphafold3' FANS OUT one AF3 job per sequence "
                        "(Aurum-only) into <gen>_af3/structs/ + af_output/, then extracts PAE "
                        "-> <gen>_af3/pae/. Both also enable global_confidence + "
                        "interdomain_pae. Ignored if --structs_dir is given.")
    p.add_argument("--af3_model_seeds", type=int, nargs="+", default=[42],
                   help="With --fold alphafold3: AF3 model seeds per sequence (default 42).")
    p.add_argument("--af3_cofold", default="none",
                   choices=["none", "mg", "mg_ppi", "mg_gpp", "mg_fpp", "mg_ggpp", "mg_gfpp", "mg_ee"],
                   help="With --fold alphafold3: co-fold the class-I TPS active site for a "
                        "HOLO prediction. 'none' (default) = apo protein only; 'mg' = the "
                        "trinuclear Mg2+ cluster; 'mg_ppi' = Mg2+ cluster + a bare diphosphate "
                        "head group; 'mg_gpp|mg_fpp|mg_ggpp|mg_gfpp' = Mg2+ cluster + ONE forced "
                        "prenyl-PP substrate (SMILES) for EVERY design; 'mg_ee' = Mg2+ cluster + "
                        "each design's EnzymeExplorer-predicted substrate (needs --enzymeexplorer_csv; "
                        "non-co-foldable EE calls fall back to Mg-only). Any non-'none' mode "
                        "ENABLES the holo tools (ion_site_check, substrate_positioning); 'none' "
                        "(or --no_holo_tools) turns co-folding AND those downstream tools off. "
                        "Filenames stay <ID>.pdb regardless.")
    p.add_argument("--enzymeexplorer_csv", default=None,
                   help="EnzymeExplorer seq-only CSV for --af3_cofold mg_ee (per-design substrate). "
                        "REQUIRED with mg_ee: the AF3 fan-out runs on the login node and cannot "
                        "afterok-wait on the in-pipeline ee_seq SLURM job, so the EE predictions "
                        "must already exist (run enzyme_explorer_sequence_only first).")
    p.add_argument("--no_holo_tools", action="store_true",
                   help="Force-skip the co-fold-dependent structure tools (ion_site_check, "
                        "substrate_positioning) even when holo structures are present (e.g. an "
                        "externally co-folded --structs_dir). By default these run iff "
                        "--af3_cofold is not 'none'.")
    p.add_argument("--no_fold_pae", action="store_true",
                   help="With --fold: do NOT save/extract the PAE matrices (skips "
                        "global_confidence + interdomain_pae). Default: save them.")
    p.add_argument("--known_structs_dir", default=None,
                   help="Known-TPS reference structures dir for the foldseek structural-identity "
                        "metric (e.g. MARTS-DB train AFDB structures, or EE reference domains).")
    p.add_argument("--pae_dir", default=None,
                   help="Dir of per-structure PAE matrices (<ID>_pae.npz from a PAE-saving "
                        "ESMFold run or the AF3 extract_pae step) -> enables global_confidence "
                        "(pTM) and interdomain_pae.")
    p.add_argument("--train_structs_dir", default=None, help="(v2) train structures dir.")
    p.add_argument("--train_embeddings_path", default=None, help="Precomputed train embeddings CSV.")
    p.add_argument("--data_colors", nargs=2, default=["dodgerblue", "goldenrod"],
                   metavar=("TRAIN", "GEN"), help="Matplotlib colors for train/generated.")
    p.add_argument("--self_consistency", action="store_true",
                   help="Back-compat alias for '--include self_consistency': add the heavy "
                        "scRMSD self-consistency step (ProteinMPNN -> ESMFold refold -> RMSD; "
                        "~1-2.5 min/structure x num_seqs on GPU). Off by default.")
    p.add_argument("--self_consistency_num_seqs", type=int, default=8,
                   help="ProteinMPNN sequences per structure for self_consistency (default 8).")
    # Single pipeline-wide top-k for the three feeders that emit *_topk.csv for the
    # knn + sdr_divergence consumers. Default comes from pipeline_tools.json
    # "_settings.top_k" (else the built-in DEFAULT_TOP_K); this overrides it.
    p.add_argument("--top_k", type=int, default=None,
                   help="Neighbours per query in the feeders' *_topk.csv (consumed by knn + "
                        f"sdr_divergence). Default: pipeline_tools.json _settings.top_k or "
                        f"{DEFAULT_TOP_K}. Only emitted when knn or sdr_divergence is enabled.")
    p.add_argument("--knn_label_file", default=DEFAULT_KNN_LABEL_FILE,
                   help="CSV mapping reference_id,label for k-NN label transfer "
                        f"(default {os.path.relpath(DEFAULT_KNN_LABEL_FILE, REPO)}).")
    p.add_argument("--knn_calibration", default=DEFAULT_KNN_CALIBRATION,
                   help="k-NN calibration JSON from `calibrate` "
                        f"(default {os.path.relpath(DEFAULT_KNN_CALIBRATION, REPO)}; "
                        "sequence+embedding+structural spaces).")
    p.add_argument("--substrate_label_file", default=DEFAULT_SUBSTRATE_LABEL_FILE,
                   help="CSV mapping reference_id,label (SUBSTRATE classes) for the "
                        "substrate_class combiner "
                        f"(default {os.path.relpath(DEFAULT_SUBSTRATE_LABEL_FILE, REPO)}).")
    p.add_argument("--substrate_calibration", default=DEFAULT_SUBSTRATE_CALIBRATION,
                   help="Substrate k-NN calibration JSON for substrate_class "
                        f"(default {os.path.relpath(DEFAULT_SUBSTRATE_CALIBRATION, REPO)}).")
    # Config-driven tool selection (precedence: config defaults -> only -> include -> exclude).
    p.add_argument("--only", default=None, metavar="A,B,...",
                   help="Run ONLY these tool keys (comma-separated); plots is kept on unless "
                        "also excluded. See --list-tools for valid keys.")
    p.add_argument("--include", default=None, metavar="A,B,...",
                   help="Force-enable these tool keys on top of the config defaults.")
    p.add_argument("--exclude", default=None, metavar="A,B,...",
                   help="Force-disable these tool keys.")
    p.add_argument("--list-tools", action="store_true", dest="list_tools",
                   help="Print each tool key (default/branch/description) and exit.")
    p.add_argument("--dry-run", action="store_true", help="Print the submission plan; don't submit.")
    args = p.parse_args()

    catalog = load_tools_catalog()
    if args.list_tools:
        width = max(len(k) for k in catalog)
        print(f"Tools catalog ({TOOLS_CONFIG if os.path.isfile(TOOLS_CONFIG) else 'built-in defaults'}):")
        print(f"  {'KEY'.ljust(width)}  DEFAULT  BRANCH      DESCRIPTION")
        for key in catalog:
            v = catalog[key]
            on = "on " if v.get("default") else "off"
            print(f"  {key.ljust(width)}  {on}      {v.get('branch','?'):<10}  {v.get('description','')}")
        return

    args.fasta_path = os.path.abspath(args.fasta_path)
    for opt in ("train_path", "structs_dir", "known_structs_dir", "pae_dir", "train_structs_dir",
                "train_embeddings_path", "knn_label_file", "knn_calibration",
                "substrate_label_file", "substrate_calibration"):
        if getattr(args, opt):
            setattr(args, opt, os.path.abspath(getattr(args, opt)))

    # AF3 folding is Aurum-only (custom partition + apptainer image). Fail fast with a
    # clear message rather than silently [cap]-skipping the producer (which would cascade
    # into the whole structure branch skipping on unmet deps).
    if args.fold == "alphafold3" and not args.structs_dir \
            and "alphafold" not in CLUSTERS[args.cluster]["caps"]:
        raise SystemExit(
            f"[FAIL] --fold alphafold3 is not available on '{args.cluster}' (AlphaFold3 is "
            "Aurum-only). Use --fold esmfold (both clusters), or pre-fold and pass --structs_dir.")

    # Resolve + export the SLURM account so every submitted sbatch inherits it via
    # $SBATCH_ACCOUNT. The per-install grant id lives in scripts/<cluster>/config.local.sh.
    account = resolve_sbatch_account(args.cluster, args.account)
    if not account and not args.dry_run:
        print(f"[warn] No SLURM account resolved for '{args.cluster}' — jobs will use your default "
              f"SLURM account, which may lack a partition association and be rejected. Set it in "
              f"scripts/{args.cluster}/config.local.sh (see config.local.sh.example), or pass "
              f"--account / export SBATCH_ACCOUNT.", file=sys.stderr)

    # Resolve the single pipeline top-k: CLI --top_k wins, else JSON _settings/DEFAULT.
    if args.top_k is None:
        args.top_k = load_top_k()

    ts = datetime.now().strftime("%Y_%m_%d_%H%M%S")
    log_dir = os.path.join(REPO, "logs", f"run_eval_pipeline-{ts}")
    if not args.dry_run:
        os.makedirs(log_dir, exist_ok=True)

    enabled = resolve_enabled_tools(catalog, args)
    steps = build_steps(args, enabled)
    print(f"=== tps_eval pipeline on '{args.cluster}' — {len(steps)} steps; logs: {log_dir} ===")
    Engine(args.cluster, log_dir, dry_run=args.dry_run).run(steps)


if __name__ == "__main__":
    main()
