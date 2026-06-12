from __future__ import annotations

import sys
from functools import reduce
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

import pandas as pd

CURRENT_DIR = Path(__file__).resolve().parent
SRC_DIR = CURRENT_DIR.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences  # noqa: E402
from data.embeddings import load_embeddings  # noqa: E402


DEFAULT_LOAD: Tuple[str, ...] = (
    "sequence",
    "max_sequence_identity",
    "max_sequence_identity_self",
    "local_sequence_search",
    "local_sequence_search_self",
    "embedding",
    "min_embedding_distance",
    "min_embedding_distance_self",
    "enzyme_explorer_sequence_only",
    "enzyme_explorer",
    "motifs",
    "soluprot",
    "motif_pair_distance",
    "esm_pseudo_perplexity",
    "swissprot_search",
)


TYPE_TO_SUBSTRATE = {
    "mono": [
        "CC(C)=CCCC(C)=CCOP([O_])(=O)OP([O_])([O_])=O",
        "CC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O (Geranyl pyrophosphate)",
    ],
    "sesq": [
        "CC(C)=CCCC(C)=CCCC(C)=CCOP([O_])(=O)OP([O_])([O_])=O",
        "CC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O (Farnesyl pyrophosphate)",
    ],
    "di": [
        "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O_])(=O)OP([O_])([O_])=O",
        "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O (Geranylgeranyl pyrophosphate)",
    ],
    "sester": [
        "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O_])(=O)OP([O_])([O_])=O",
        "CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP([O-])(=O)OP([O-])([O-])=O (Geranylfarnesyl pyrophosphate)",
    ],
    "tri": [
        "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC1OC1(C)C",
        "CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC1OC1(C)C ((S)-2,3-epoxysqualene)",
    ],
}


def get_target_substrate(tps_type: str, df_names: Iterable[str]) -> Optional[str]:
    target_substrate_names = TYPE_TO_SUBSTRATE[tps_type]
    df_names_set = set(df_names)
    for name in target_substrate_names:
        if name in df_names_set:
            return name
    return None


def construct_result_paths(fasta_path: str) -> Tuple[str, ...]:
    if fasta_path.endswith(".fasta"):
        remove_length = 6
    elif fasta_path.endswith(".fa"):
        remove_length = 3
    else:
        raise ValueError("The fasta file path must end with '.fasta' or '.fa'.")
    partial = fasta_path[: -remove_length]

    return (
        fasta_path,
        partial + "_max_sequence_identity.csv",
        partial + "_max_sequence_identity_self.csv",
        partial + "_embedding_esm1b.csv",
        partial + "_embedding_esm1b_min_embedding_distance.csv",
        partial + "_embedding_esm1b_min_embedding_distance_self.csv",
        partial + "_enzyme_explorer_sequence_only.csv",
        partial + "_enzyme_explorer.csv",
        partial + "_motifs.csv",
        partial + "_soluprot.csv",
    )


def _fasta_partial(fasta_path: str) -> str:
    """The `<input>` stem (path without the .fasta/.fa extension) that metric
    tools prepend to their CSV save names (e.g. `<input>_motif_pair_distance.csv`)."""
    if fasta_path.endswith(".fasta"):
        return fasta_path[:-6]
    if fasta_path.endswith(".fa"):
        return fasta_path[:-3]
    raise ValueError("The fasta file path must end with '.fasta' or '.fa'.")


def _strip_column_names(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={col: col.strip() for col in df.columns})


def _load_enzyme_explorer(
    csv_path: str,
    *,
    tps_type: Optional[str],
    rename_isTPS_to: Optional[str] = None,
) -> pd.DataFrame:
    scores = pd.read_csv(csv_path)
    scores = _strip_column_names(scores)

    # --- Schema normalization (old EnzymeExplorer vs the `revision` branch) ---
    # Old EE: "ID" + an "isTPS" probability column.
    # Revision EE (predict_sequences_only / predict_with_structures): lowercase
    # "id" + per-class "<class>_score" / "<class>_p_calibrated"; the TPS-vs-not
    # probability is "TPS_p_calibrated" (calibrated, preferred) or "TPS_score".
    if "ID" not in scores.columns and "id" in scores.columns:
        scores = scores.rename(columns={"id": "ID"})
    if "isTPS" not in scores.columns:
        if "TPS_p_calibrated" in scores.columns:
            scores = scores.rename(columns={"TPS_p_calibrated": "isTPS"})
        elif "TPS_score" in scores.columns:
            scores = scores.rename(columns={"TPS_score": "isTPS"})
        else:
            raise KeyError(
                f"{csv_path}: no TPS-probability column found "
                "(expected 'isTPS', 'TPS_p_calibrated', or 'TPS_score')"
            )

    columns = ["ID", "isTPS"]
    target_substrate: Optional[str] = None
    if tps_type is not None:
        target_substrate = get_target_substrate(tps_type, scores.columns)
        if target_substrate is not None:
            columns.append(target_substrate)

    scores = scores.loc[:, columns].copy()
    if target_substrate is not None:
        scores = scores.rename(columns={target_substrate: "target_substrate"})
    if rename_isTPS_to is not None:
        scores = scores.rename(columns={"isTPS": rename_isTPS_to})
    return scores


def _strip_id_column(df: pd.DataFrame) -> pd.DataFrame:
    if "ID" in df.columns:
        print("Warning: stripping ID column of parts after whitespaces.")
        df = df.copy()
        df["ID"] = df["ID"].astype(str).map(lambda x: x.split(" ", 1)[0])
    return df


def _outer_join_on_id(frames: List[pd.DataFrame]) -> pd.DataFrame:
    if not frames:
        # No input CSVs survived selection (e.g. a structure-only run with no
        # sequence-branch outputs). Return an empty frame that still carries the
        # "ID" key so downstream merge/column access degrades gracefully instead
        # of crashing reduce() with an empty iterable.
        return pd.DataFrame(columns=["ID"])
    if len(frames) == 1:
        return frames[0]
    return reduce(lambda left, right: left.merge(right, on="ID", how="outer"), frames)


def load_results(
    fasta_path: str,
    max_sequence_identity_path: Optional[str] = None,
    max_sequence_identity_self_path: Optional[str] = None,
    embedding_esm1b_path: Optional[str] = None,
    min_embedding_distance_esm1b_path: Optional[str] = None,
    min_embedding_distance_esm1b_self_path: Optional[str] = None,
    enzyme_explorer_sequence_only_path: Optional[str] = None,
    enzyme_explorer_path: Optional[str] = None,
    motifs_path: Optional[str] = None,
    soluprot_path: Optional[str] = None,
    *,
    tps_type: Optional[str] = None,
    load: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Load combined results for a fasta file.

    If only `fasta_path` is provided, all auxiliary CSV paths are derived from it
    using `construct_result_paths`.
    """
    if any(
        p is None
        for p in (
            max_sequence_identity_path,
            max_sequence_identity_self_path,
            embedding_esm1b_path,
            min_embedding_distance_esm1b_path,
            min_embedding_distance_esm1b_self_path,
            enzyme_explorer_sequence_only_path,
            enzyme_explorer_path,
            motifs_path,
            soluprot_path,
        )
    ):
        (
            fasta_path,
            max_sequence_identity_path,
            max_sequence_identity_self_path,
            embedding_esm1b_path,
            min_embedding_distance_esm1b_path,
            min_embedding_distance_esm1b_self_path,
            enzyme_explorer_sequence_only_path,
            enzyme_explorer_path,
            motifs_path,
            soluprot_path,
        ) = construct_result_paths(fasta_path)

    if load is None:
        load = DEFAULT_LOAD
    load_set = set(load)

    frames: List[pd.DataFrame] = []

    if "sequence" in load_set:
        records = load_fasta_sequences(fasta_path, remove_padding=True, load_identifiers=True)
        ids = [rid for rid, _ in records]
        seqs = [seq for _, seq in records]
        msa = load_fasta_sequences(fasta_path, remove_padding=False)
        sequences_df = pd.DataFrame({"ID": ids, "sequence": seqs, "msa": msa})
        frames.append(sequences_df)

    if "embedding" in load_set:
        embedding_df = load_embeddings(embedding_esm1b_path)
        frames.append(embedding_df)

    if "enzyme_explorer_sequence_only" in load_set:
        df = _load_enzyme_explorer(
            enzyme_explorer_sequence_only_path,
            tps_type=tps_type,
            rename_isTPS_to="isTPS_seq",
        )
        frames.append(df)

    if "enzyme_explorer" in load_set:
        df = _load_enzyme_explorer(
            enzyme_explorer_path,
            tps_type=tps_type,
        )
        frames.append(df)

    if "max_sequence_identity" in load_set:
        df = pd.read_csv(max_sequence_identity_path)
        frames.append(df)

    if "max_sequence_identity_self" in load_set:
        df = pd.read_csv(max_sequence_identity_self_path)
        df = df.rename(
            columns={
                "sequence_identity": "sequence_identity_self",
                "sequence_similarity": "sequence_similarity_self",
            }
        )
        frames.append(df)

    if "min_embedding_distance" in load_set:
        df = pd.read_csv(min_embedding_distance_esm1b_path)
        frames.append(df)

    if "min_embedding_distance_self" in load_set:
        df = pd.read_csv(min_embedding_distance_esm1b_self_path)
        df = df.rename(
            columns={"min_embedding_distance": "min_embedding_distance_self"}
        )
        frames.append(df)

    if "motifs" in load_set:
        df = pd.read_csv(motifs_path)
        if "sequence" in df.columns:
            df = df.drop(columns=["sequence"])
        frames.append(df)

    if "soluprot" in load_set:
        df = pd.read_csv(soluprot_path)
        df = df.rename(columns={"fa_id": "ID"})
        frames.append(df)

    # --- Newer sequence-branch metric CSVs (keyed by ID, sibling of the fasta) ---
    # These follow the standard `<input>_<tool>.csv` naming, so derive the path
    # straight from the fasta stem rather than threading more positional args.
    if "motif_pair_distance" in load_set:
        df = pd.read_csv(_fasta_partial(fasta_path) + "_motif_pair_distance.csv")
        df = _strip_column_names(df)
        # Drop the raw motif-string columns; keep only the scalar metrics + ID.
        keep = [c for c in df.columns if c not in ("ddxxd_motif", "nse_dte_motif")]
        frames.append(df.loc[:, keep])

    if "esm_pseudo_perplexity" in load_set:
        df = pd.read_csv(_fasta_partial(fasta_path) + "_esm_pseudo_perplexity.csv")
        df = _strip_column_names(df)
        if "n_residues" in df.columns:
            df = df.drop(columns=["n_residues"])
        frames.append(df)

    if "swissprot_search" in load_set:
        df = pd.read_csv(_fasta_partial(fasta_path) + "_swissprot_search.csv")
        df = _strip_column_names(df)
        # `swissprot_top_hit` is the accession string — not a plottable metric.
        if "swissprot_top_hit" in df.columns:
            df = df.drop(columns=["swissprot_top_hit"])
        frames.append(df)

    # Fast LOCAL (MMseqs2) sequence search. gen-vs-train -> default file; within-set
    # (self) -> the _self file written by the pipeline's local_search_{tag} step, with
    # columns renamed to *_self so the two frames don't collide on the outer join.
    if "local_sequence_search" in load_set:
        df = pd.read_csv(_fasta_partial(fasta_path) + "_local_sequence_search.csv")
        df = _strip_column_names(df)
        frames.append(df)

    if "local_sequence_search_self" in load_set:
        df = pd.read_csv(_fasta_partial(fasta_path) + "_local_sequence_search_self.csv")
        df = _strip_column_names(df)
        df = df.rename(
            columns={
                "local_sequence_identity": "local_sequence_identity_self",
                "local_sequence_similarity": "local_sequence_similarity_self",
                "local_coverage": "local_coverage_self",
            }
        )
        frames.append(df)

    # k-NN coarse-label transfer (gen-only). `confidence` (numeric) + `predicted_label`
    # (categorical) are the plottable columns; the rest (per-space votes) are diagnostics.
    if "knn_label_transfer" in load_set:
        df = pd.read_csv(_fasta_partial(fasta_path) + "_knn_label_transfer.csv")
        df = _strip_column_names(df)
        frames.append(df)

    # Substrate-class combiner (gen-only). Rename the final `confidence` to
    # `substrate_confidence` so it doesn't collide with the k-NN `confidence` target.
    if "substrate_class" in load_set:
        df = pd.read_csv(_fasta_partial(fasta_path) + "_substrate_class.csv")
        df = _strip_column_names(df)
        if "confidence" in df.columns:
            df = df.rename(columns={"confidence": "substrate_confidence"})
        frames.append(df)

    frames = [_strip_id_column(f) for f in frames]

    return _outer_join_on_id(frames)


def load_structure_results(csv_path: str) -> pd.DataFrame:
    """Load a single structure-metric CSV (keyed by ID).

    Structure metrics are saved as `<structs_dir>_<tool>.csv` next to the
    structures directory — i.e. they are NOT derivable from a fasta stem the
    way the sequence metrics are. The plot layer discovers these files by name
    in the input directory and loads them through this helper. The bookkeeping
    `n_residues` column (present in several of them) is dropped so it isn't
    mistaken for a metric. Raises FileNotFoundError when the CSV is absent, so
    the caller's skip-missing-input behavior fires as for any other target.
    """
    df = pd.read_csv(csv_path)
    df = _strip_column_names(df)
    if "n_residues" in df.columns:
        df = df.drop(columns=["n_residues"])
    return _strip_id_column(df)
