from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple

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


def _extract_pae(output, seq_len: int):
    """Pull the (L, L) Predicted Aligned Error matrix (Angstrom) out of an
    EsmForProteinFolding output, or return None if the head is unavailable.

    HF's output carries ``predicted_aligned_error`` as (batch, L, L) in Angstrom
    (the bin-expectation over ``aligned_confidence_probs``). We take batch element
    0 and crop to the real sequence length defensively (ESMFold here folds with
    add_special_tokens=False, so no BOS/EOS padding, but cropping is harmless and
    robust to a future change). Returns a float32 numpy (L, L) or None — saving PAE
    must DEGRADE GRACEFULLY if the field is missing (the .pdb is unaffected)."""
    import numpy as np

    pae = None
    if isinstance(output, dict):
        pae = output.get("predicted_aligned_error")
    else:
        pae = getattr(output, "predicted_aligned_error", None)
    if pae is None:
        return None
    arr = pae.detach().to("cpu").float().numpy()
    if arr.ndim == 3:  # (batch, L, L) -> first element
        arr = arr[0]
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        return None
    # Crop to the actual residue count if any special tokens slipped in.
    L = min(seq_len, arr.shape[0])
    return np.ascontiguousarray(arr[:L, :L], dtype=np.float32)


def _extract_ptm(output):
    """Pull the global pTM (predicted TM-score; 0-1, higher=better — whole-fold
    confidence) scalar out of an EsmForProteinFolding output, or None if absent.

    HF's output carries ``ptm`` as a 0-d (or batch) tensor. We take the scalar (batch
    element 0 if batched) and return a python float — or None so saving DEGRADES
    GRACEFULLY (the .pdb and the PAE fields are unaffected). ESMFold is single-chain,
    so there is no ipTM."""
    ptm = output.get("ptm") if isinstance(output, dict) else getattr(output, "ptm", None)
    if ptm is None:
        return None
    try:
        arr = ptm.detach().to("cpu").float().numpy()
        return float(arr.reshape(-1)[0])
    except Exception:  # noqa: BLE001 - any odd shape/type -> treat as unavailable
        return None


def fold_sequence(
    model, tokenizer, device: str, sequence: str, *, chunk_size: Optional[int] = None
) -> Tuple[str, "object", "object"]:
    """Fold a single sequence. Returns ``(pdb_str, pae, ptm)`` where ``pdb_str`` has
    the per-residue pLDDT in the B-factor field rescaled to the AlphaFold 0-100
    convention, ``pae`` is the (L, L) float32 Predicted-Aligned-Error matrix in
    Angstrom (or None if the PAE head is unavailable), and ``ptm`` is the global pTM
    fold-confidence scalar (float in 0-1, or None if unavailable). The PDB string is
    BYTE-FOR-BYTE identical to the pre-PAE behaviour — PAE/pTM extraction is
    read-only."""
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
    pdb_str = _rescale_bfactor_to_0_100(model.output_to_pdb(output)[0])
    pae = _extract_pae(output, len(sequence))
    ptm = _extract_ptm(output)
    return pdb_str, pae, ptm


def _default_pae_dir(save_dir: str) -> str:
    """Sibling PAE dir next to the structs dir: ``<save_dir>_pae/``."""
    d = save_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_pae")


def _residue_ids_from_pdb(pdb_str: str) -> List[int]:
    """The PDB author residue number of every CA-bearing standard residue, in order
    — this is the PAE matrix's residue axis. We read it back from the PDB we just
    wrote (rather than assuming 1..L) so the npz ``residue_ids`` axis exactly matches
    the structure's numbering, the same numbering EnzymeExplorer's domain detector
    reports via PyMOL ``resi`` — keeping the consumer numbering-robust."""
    res_ids: List[int] = []
    seen = set()
    for line in pdb_str.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            key = (line[21], resseq, line[26])  # chain, resSeq, iCode
            if key not in seen:
                seen.add(key)
                res_ids.append(resseq)
    return res_ids


def _save_pae(pae, pdb_str: str, pae_path: str, *, seq_len: int, ptm=None) -> bool:
    """Write the shared ``<ID>_pae.npz`` (schema: pae, residue_ids, n_residues,
    source, ptm). Returns True if written, False if `pae` was None (head
    unavailable). The schema is IDENTICAL to the AF3 extractor so interdomain_pae and
    global_confidence are source-agnostic. residue_ids is read from the PDB CA
    records; if that count disagrees with the matrix it falls back to a 1..L axis
    (ESMFold numbers contiguously). ``ptm`` is the global fold-confidence scalar
    (0-1), stored as a float32 (NaN if unavailable — additive, never affects the
    existing pae/residue_ids/n_residues/source fields). ESMFold is single-chain, so
    no ipTM is stored."""
    import numpy as np

    if pae is None:
        return False
    res_ids = _residue_ids_from_pdb(pdb_str)
    if len(res_ids) != pae.shape[0]:
        res_ids = list(range(1, pae.shape[0] + 1))
    os.makedirs(os.path.dirname(os.path.abspath(pae_path)), exist_ok=True)
    np.savez_compressed(
        pae_path,
        pae=np.ascontiguousarray(pae, dtype=np.float32),
        residue_ids=np.asarray(res_ids, dtype=np.int32),
        n_residues=np.int64(pae.shape[0]),
        source="esmfold",
        ptm=np.float32(np.nan if ptm is None else ptm),
    )
    return True


def fold_fasta(
    fasta_path: str,
    save_dir: str,
    *,
    skip_existing: bool = True,
    chunk_size: Optional[int] = None,
    device: Optional[str] = None,
    save_pae: bool = True,
    pae_dir: Optional[str] = None,
) -> List[str]:
    """Fold every sequence in `fasta_path` with ESMFold, writing `<ID>.pdb` into
    `save_dir`. Returns the list of written PDB paths. The output dir mirrors the
    AlphaFold structs/ layout and is consumed unchanged by run_plddt /
    run_structural_identity (ID = FASTA record id = filename stem).

    When `save_pae` (default True), ALSO writes the (L, L) Predicted-Aligned-Error
    matrix per sequence as ``<ID>_pae.npz`` into `pae_dir` (default
    ``<save_dir>_pae/``), in the shared schema consumed by
    structure_metrics/interdomain_pae. PAE saving is ADDITIVE — the ``<ID>.pdb`` is
    unchanged — and degrades gracefully (warns, keeps the pdb) if the PAE head is
    unavailable for a sequence."""
    identifiers, sequences = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    os.makedirs(save_dir, exist_ok=True)
    if save_pae and pae_dir is None:
        pae_dir = _default_pae_dir(save_dir)
    if save_pae:
        os.makedirs(pae_dir, exist_ok=True)

    # Decide which records still need folding before paying the model-load cost.
    # A record is "done" only when its .pdb exists AND (if requested) its PAE npz
    # exists — so enabling PAE on a previously-folded structs dir re-folds to fill
    # in the missing matrices.
    todo = []
    for record_id, sequence in zip(identifiers, sequences):
        stem = _sanitize_id(record_id)
        out_path = os.path.join(save_dir, stem + ".pdb")
        pae_path = os.path.join(pae_dir, stem + "_pae.npz") if save_pae else None
        if skip_existing and os.path.isfile(out_path) and (
            not save_pae or os.path.isfile(pae_path)
        ):
            print(f"[skip] {stem}: {out_path} already exists")
            continue
        todo.append((stem, sequence, out_path, pae_path))

    if not todo:
        print("Nothing to fold (all outputs already exist).")
        return []

    tokenizer, model, device = _load_model(device)

    written: List[str] = []
    n = len(todo)
    n_pae = 0
    for i, (stem, sequence, out_path, pae_path) in enumerate(todo, start=1):
        print(f"[{i}/{n}] folding {stem} ({len(sequence)} aa) -> {out_path}")
        try:
            pdb_str, pae, ptm = fold_sequence(model, tokenizer, device, sequence, chunk_size=chunk_size)
        except Exception as exc:  # noqa: BLE001 - keep folding the rest of the batch
            print(f"  [warn] failed to fold {stem}: {exc}")
            continue
        with open(out_path, "w") as fh:
            fh.write(pdb_str)
        written.append(out_path)
        if save_pae:
            if _save_pae(pae, pdb_str, pae_path, seq_len=len(sequence), ptm=ptm):
                n_pae += 1
            else:
                print(f"  [warn] {stem}: PAE head unavailable; wrote .pdb only (no {os.path.basename(pae_path)})")

    print(f"Wrote {len(written)}/{n} structure(s) to {save_dir}")
    if save_pae:
        print(f"Wrote {n_pae}/{len(written)} PAE matrix file(s) to {pae_dir}")
    return written
