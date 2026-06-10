from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import List, Optional

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers

# ESMFold (facebook/esmfold_v1, via HuggingFace transformers) folds a single-chain
# protein sequence into a PDB string whose B-factor field holds the per-residue
# pLDDT. HF's output_to_pdb writes that pLDDT on a 0-1 scale; we rescale it to
# 0-100 (AlphaFold convention) on write so the output is drop-in. We write one <ID>.pdb per FASTA
# record into a structs dir, mirroring the AlphaFold structs/ layout, so the SAME
# downstream tools (structure_metrics/plddt, foldseek/structural_identity) consume
# the output unchanged (ID = FASTA record id = filename stem).
MODEL_NAME = "facebook/esmfold_v1"

# Sequences longer than this trigger chunked attention to keep peak GPU memory
# bounded. ESMFold's folding trunk has O(L^2) attention; chunking trades a little
# speed for a lot of memory headroom on long sequences. TPS sequences are usually
# 300-600 aa (fit comfortably on a single A100 without chunking).
CHUNK_LENGTH_THRESHOLD = 600
DEFAULT_CHUNK_SIZE = 64


def _sanitize_id(record_id: str) -> str:
    """FASTA ids can carry path separators / whitespace that break a filename.
    Keep only the leading whitespace-stripped token's basename-safe form."""
    stem = record_id.strip().split()[0] if record_id.strip() else record_id
    return stem.replace(os.sep, "_")


def _load_model(device: Optional[str] = None):
    """Load the ESMFold model onto the best available device. Imported lazily so
    that importing this module (e.g. for tests) doesn't pull in torch/transformers."""
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Loading {MODEL_NAME} onto {device} ...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForProteinFolding.from_pretrained(MODEL_NAME)
    model = model.to(device)
    model.eval()
    # esm_s embeddings can be computed in float16 on GPU to halve the backbone's
    # memory; the folding trunk stays float32 for numerical stability.
    if device == "cuda":
        model.esm = model.esm.half()
    return tokenizer, model, device


def _rescale_bfactor_to_0_100(pdb_str: str) -> str:
    """ESMFold's HF `output_to_pdb` writes pLDDT into the B-factor field on a 0-1
    scale (`categorical_lddt` returns the bin-expectation in [0,1], passed straight
    into `b_factors` with no x100). AlphaFold — and tps_eval's downstream
    structure_metrics/plddt (CONFIDENT_THRESHOLD=70) — expect pLDDT on a 0-100
    scale. Rescale the B-factor column (PDB cols 61-66, 0-based 60:66) by x100 so
    the output is drop-in compatible with the AlphaFold structs/ layout and the
    existing pLDDT tool consumes it unchanged."""
    out_lines = []
    for line in pdb_str.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                b = float(line[60:66]) * 100.0
                line = line[:60] + f"{b:6.2f}" + line[66:]
            except ValueError:
                pass  # leave malformed/blank B-factor columns untouched
        out_lines.append(line)
    return "\n".join(out_lines) + "\n"


def fold_sequence(model, tokenizer, device: str, sequence: str, *, chunk_size: Optional[int] = None) -> str:
    """Fold a single sequence into a PDB string (pLDDT in the B-factor field,
    rescaled to the AlphaFold 0-100 convention)."""
    import torch

    if chunk_size is not None:
        model.trunk.set_chunk_size(chunk_size)
    elif len(sequence) > CHUNK_LENGTH_THRESHOLD:
        model.trunk.set_chunk_size(DEFAULT_CHUNK_SIZE)
    else:
        model.trunk.set_chunk_size(None)

    tokenized = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
    tokenized = {k: v.to(device) for k, v in tokenized.items()}
    with torch.no_grad():
        output = model(**tokenized)
    # output_to_pdb returns one PDB string per batch element with the per-residue
    # pLDDT in the B-factor field on a 0-1 scale; rescale to 0-100 (AlphaFold
    # convention) so the structs dir flows through structure_metrics/plddt unchanged.
    pdb_str = model.output_to_pdb(output)[0]
    return _rescale_bfactor_to_0_100(pdb_str)


def fold_fasta(
    fasta_path: str,
    save_dir: str,
    *,
    skip_existing: bool = True,
    chunk_size: Optional[int] = None,
    device: Optional[str] = None,
) -> List[str]:
    """Fold every sequence in `fasta_path` with ESMFold, writing `<ID>.pdb` into
    `save_dir`. Returns the list of written PDB paths. The output dir mirrors the
    AlphaFold structs/ layout and is consumed unchanged by run_plddt /
    run_structural_identity (ID = FASTA record id = filename stem)."""
    identifiers, sequences = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    os.makedirs(save_dir, exist_ok=True)

    # Decide which records still need folding before paying the model-load cost.
    todo = []
    for record_id, sequence in zip(identifiers, sequences):
        stem = _sanitize_id(record_id)
        out_path = os.path.join(save_dir, stem + ".pdb")
        if skip_existing and os.path.isfile(out_path):
            print(f"[skip] {stem}: {out_path} already exists")
            continue
        todo.append((stem, sequence, out_path))

    if not todo:
        print("Nothing to fold (all outputs already exist).")
        return []

    tokenizer, model, device = _load_model(device)

    written: List[str] = []
    n = len(todo)
    for i, (stem, sequence, out_path) in enumerate(todo, start=1):
        print(f"[{i}/{n}] folding {stem} ({len(sequence)} aa) -> {out_path}")
        try:
            pdb_str = fold_sequence(model, tokenizer, device, sequence, chunk_size=chunk_size)
        except Exception as exc:  # noqa: BLE001 - keep folding the rest of the batch
            print(f"  [warn] failed to fold {stem}: {exc}")
            continue
        with open(out_path, "w") as fh:
            fh.write(pdb_str)
        written.append(out_path)

    print(f"Wrote {len(written)}/{n} structure(s) to {save_dir}")
    return written
