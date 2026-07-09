#!/usr/bin/env python3 -u
#
# Matous Soldat, 2026
#
# Compute SaProt (structure-aware PLM) per-protein embeddings.
#
# SaProt consumes a "structure-aware" (SA) token sequence: per residue, the amino
# acid letter is concatenated with its foldseek 3Di structural token (e.g. "Mp" =
# residue M in 3Di-state p). We build that SA sequence for each input PDB structure
# via SaProt's own get_struc_seq utility (which shells out to foldseek), run SaProt,
# mean-pool the per-residue embeddings, and emit one vector per structure keyed by
# the PDB filename stem (== Enzyme_marts_ID for the MARTS-DB structures).
#
# Mirrors the conventions of src/esm/extract_embeddings.py: first column "id", then
# embedding dims 0..D-1; mean-pooled over real residues (excluding BOS/EOS).
#
# Model:    westlake-repl/SaProt_650M_AF2 (ESM-2 650M backbone, vocab = AA x 3Di)
# foldseek: produced by SaProt's foldseek_util.get_struc_seq

import os
import sys
import glob
import time
import argparse

import torch
import pandas as pd

from transformers import EsmTokenizer, EsmForMaskedLM


DEFAULT_MODEL = "westlake-repl/SaProt_650M_AF2"


def create_parser():
    parser = argparse.ArgumentParser(
        description="Extract mean-pooled SaProt structure-aware per-protein embeddings."
    )
    parser.add_argument(
        "--structs_dir",
        type=str,
        required=True,
        help="Directory of PDB structures named <ID>.pdb (ID becomes the row key).",
    )
    parser.add_argument(
        "--foldseek",
        type=str,
        required=True,
        help="Path to the foldseek binary (used to derive the 3Di structural tokens).",
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        required=True,
        help="Output CSV path. First column 'id', then embedding dims 0..D-1.",
    )
    parser.add_argument(
        "--saprot_repo",
        type=str,
        default=None,
        help="Path to a cloned SaProt repo (for utils/foldseek_util.get_struc_seq). "
        "If omitted, expects 'utils.foldseek_util' to be importable on sys.path.",
    )
    parser.add_argument(
        "--ids_csv",
        type=str,
        default=None,
        help="Optional CSV with an 'Enzyme_marts_ID' column. If given, only structures "
        "whose stem is in that column are embedded (and coverage is reported against it).",
    )
    parser.add_argument(
        "--id_column",
        type=str,
        default="Enzyme_marts_ID",
        help="Column name in --ids_csv holding the structure stems / row keys.",
    )
    parser.add_argument(
        "--model_location",
        type=str,
        default=DEFAULT_MODEL,
        help="HuggingFace model id or local path for the SaProt model.",
    )
    parser.add_argument(
        "--chain",
        type=str,
        default="A",
        help="Chain id to extract from each PDB (foldseek key). Default 'A'.",
    )
    parser.add_argument(
        "--truncation_seq_length",
        type=int,
        default=1024,
        help="Truncate SA sequences longer than this many residues.",
    )
    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available.")
    return parser


def load_get_struc_seq(saprot_repo):
    """Import SaProt's get_struc_seq (shells out to foldseek to produce 3Di tokens)."""
    if saprot_repo is not None:
        sys.path.insert(0, saprot_repo)
    from utils.foldseek_util import get_struc_seq  # noqa: E402

    return get_struc_seq


def build_sa_sequence(get_struc_seq, foldseek, pdb_path, chain):
    """Return the combined SA sequence (AA + 3Di interleaved) for one structure.

    plddt_mask="auto": SaProt masks low-pLDDT residues only when the PDB self-identifies
    as AlphaFold output. ESMFold PDBs do not, so no masking is applied here.
    """
    seq_dict = get_struc_seq(foldseek, pdb_path, chains=[chain], plddt_mask="auto")
    if chain in seq_dict:
        _, _, combined = seq_dict[chain]
    else:
        # Fall back to the first (only) chain foldseek reported.
        first = next(iter(seq_dict.values()))
        _, _, combined = first
    return combined


def run(args):
    get_struc_seq = load_get_struc_seq(args.saprot_repo)

    print(f"Loading SaProt model: {args.model_location}", flush=True)
    tokenizer = EsmTokenizer.from_pretrained(args.model_location)
    model = EsmForMaskedLM.from_pretrained(args.model_location)
    model.eval()
    use_gpu = torch.cuda.is_available() and not args.nogpu
    if use_gpu:
        model = model.cuda()
        print("Transferred model to GPU", flush=True)

    # Collect target IDs.
    wanted = None
    if args.ids_csv is not None:
        ref = pd.read_csv(args.ids_csv)
        wanted = set(ref[args.id_column].astype(str).unique())
        print(f"Reference set: {len(wanted)} unique {args.id_column}", flush=True)

    pdb_files = sorted(glob.glob(os.path.join(args.structs_dir, "*.pdb")))
    print(f"Found {len(pdb_files)} PDB files in {args.structs_dir}", flush=True)

    embedding_labels = []
    embedding_rows = []
    failed = []

    with torch.no_grad():
        for n, pdb_path in enumerate(pdb_files):
            stem = os.path.splitext(os.path.basename(pdb_path))[0]
            if wanted is not None and stem not in wanted:
                continue
            try:
                t0 = time.time()
                sa_seq = build_sa_sequence(get_struc_seq, args.foldseek, pdb_path, args.chain)
                if not sa_seq:
                    raise ValueError("empty SA sequence")
                # SA seq is 2 chars per residue; truncate by residue count.
                max_chars = 2 * args.truncation_seq_length
                if len(sa_seq) > max_chars:
                    sa_seq = sa_seq[:max_chars]

                inputs = tokenizer(sa_seq, return_tensors="pt")
                if use_gpu:
                    inputs = {k: v.cuda() for k, v in inputs.items()}

                out = model.esm(**inputs)
                hidden = out.last_hidden_state  # (1, L+2, D)
                # Mean-pool over real residues, excluding BOS (idx 0) and EOS (idx -1).
                emb = hidden[0, 1:-1, :].mean(0).to("cpu").float().numpy()

                embedding_labels.append(stem)
                embedding_rows.append(emb)
                if (len(embedding_labels)) % 50 == 0:
                    print(
                        f"[{len(embedding_labels)}] {stem} "
                        f"(L={len(sa_seq)//2}, {time.time()-t0:.2f}s)",
                        flush=True,
                    )
            except Exception as e:  # noqa: BLE001
                failed.append((stem, str(e)))
                print(f"FAILED {stem}: {e}", flush=True)

    if not embedding_rows:
        raise RuntimeError("No embeddings were produced.")

    label_df = pd.DataFrame({"id": embedding_labels})
    embedding_df = pd.DataFrame(embedding_rows)
    embedding_df.columns = [str(c) for c in embedding_df.columns]
    df = pd.concat([label_df, embedding_df], axis=1)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_csv)), exist_ok=True)
    df.to_csv(args.output_csv, index=False)

    dim = embedding_df.shape[1]
    print(f"\nWrote {len(df)} embeddings (dim={dim}) to {args.output_csv}", flush=True)
    if wanted is not None:
        covered = len(set(embedding_labels) & wanted)
        print(
            f"Coverage: {covered}/{len(wanted)} reference {args.id_column} embedded.",
            flush=True,
        )
    if failed:
        print(f"{len(failed)} structures failed:", flush=True)
        for stem, err in failed[:20]:
            print(f"  {stem}: {err}", flush=True)


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
