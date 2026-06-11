"""Extract the PAE matrix from an AlphaFold3 af_output tree into the shared
``<ID>_pae.npz`` schema consumed by src/structure_metrics/interdomain_pae.py.

AF3 writes, per job, ``af_output/<job>/<job>_confidences.json`` whose ``"pae"``
field is the token x token Predicted Aligned Error matrix (Angstrom; PAE[i, j] is
the expected position error of token j when the prediction is aligned on token i).
The companion ``"token_res_ids"`` array gives the residue number of each token and
``"token_chain_ids"`` the chain — for a single-chain protein this is one token per
residue, so ``token_res_ids`` IS the PAE residue axis (matching ESMFold's
``<ID>_pae.npz``). ID = the job subfolder name (== the structs/ stem the rest of the
pipeline keys off, mirroring run_plddt's af3 layout detection).

AF3 also emits the global ``ptm`` (predicted TM-score; whole-fold confidence, 0-1,
higher=better) and ``iptm`` (interface pTM, multi-chain only) scalars. These live in
the per-job ``<job>_summary_confidences.json`` (the summary file) — the large
``<job>_confidences.json`` carries the matrices (``pae`` etc.) but not the scalars,
so we read the summary file alongside it when present. ``iptm`` is null/absent for
single-chain TPS designs -> stored as NaN.

Shared npz schema (identical for ESMFold and AF3 so the consumer is source-agnostic):
  * ``pae``         : float32 (L, L), Angstrom.
  * ``residue_ids`` : int32 (L,), the PDB author residue number of each PAE row/col.
  * ``n_residues``  : int scalar (== L).
  * ``source``      : str, here 'alphafold3'.
  * ``ptm``         : float32 scalar, global fold confidence (0-1; NaN if absent).
  * ``iptm``        : float32 scalar, interface pTM (multi-chain; NaN if absent).

Usage:
    python extract_pae.py --af_output <af_output dir> --pae_dir <out dir>
or point --structs_dir at the structs/ dir and we look for its sibling af_output/.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np


def _resolve_af_output(af_output: str | None, structs_dir: str | None) -> Path:
    """Locate the af_output tree from --af_output or a --structs_dir sibling.

    Mirrors run_plddt's af3 convention: an af_output dir holds one subfolder per
    job, each with ``<job>/<job>_confidences.json``. If --structs_dir is given we
    accept either the af_output dir itself or a dir whose sibling/child is
    ``af_output/``."""
    if af_output:
        return Path(af_output)
    if not structs_dir:
        raise ValueError("Pass either --af_output or --structs_dir.")
    sd = Path(structs_dir)
    for cand in (sd, sd / "af_output", sd.parent / "af_output"):
        if cand.is_dir() and any(
            (cand / e / f"{e}_confidences.json").is_file() for e in os.listdir(cand)
        ):
            return cand
    raise FileNotFoundError(
        f"No af_output tree (with <job>/<job>_confidences.json) found from {structs_dir}"
    )


def _coerce_scalar(value) -> float:
    """AF3 emits ptm/iptm as a number (or null). iptm can also be a per-chain list in
    some outputs — collapse a list to its mean. Returns float NaN when absent."""
    if value is None:
        return float("nan")
    if isinstance(value, (list, tuple)):
        nums = [float(v) for v in value if v is not None]
        return float(np.mean(nums)) if nums else float("nan")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _load_ptm_iptm(confidences_json: str, job_id: str):
    """Read the global ``ptm`` and ``iptm`` scalars for a job.

    The scalars live in the sibling ``<job>_summary_confidences.json`` (the summary
    file), not the big matrix ``<job>_confidences.json``. We also fall back to the
    matrix file itself in case a given AF3 version inlined them there. Returns
    ``(ptm, iptm)`` as floats (NaN when absent — e.g. iptm for single-chain TPS).
    """
    p = Path(confidences_json)
    summary = p.with_name(p.name.replace("_confidences.json", "_summary_confidences.json"))
    for src in (summary, p):
        if not src.is_file():
            continue
        try:
            with open(src) as fh:
                payload = json.load(fh)
        except Exception:  # noqa: BLE001 - unreadable summary -> try next / NaN
            continue
        if "ptm" in payload or "iptm" in payload:
            return _coerce_scalar(payload.get("ptm")), _coerce_scalar(payload.get("iptm"))
    print(f"  [warn] {job_id}: no ptm/iptm found (summary_confidences.json missing field) -> NaN")
    return float("nan"), float("nan")


def confidences_to_npz(confidences_json: str, out_npz: str, *, job_id: str) -> int:
    """Read one AF3 ``<job>_confidences.json`` and write the shared ``<ID>_pae.npz``.

    Returns L (the matrix dimension). For multi-chain complexes ``token_res_ids``
    can repeat across chains; we keep the raw per-token residue ids as the axis
    (the consumer maps EE per-structure PDB resi onto them) and warn — TPS designs
    here are single-chain, where token_res_ids is the residue axis 1..L."""
    with open(confidences_json) as fh:
        payload = json.load(fh)
    if "pae" not in payload:
        raise KeyError(f"No 'pae' field in {confidences_json}")
    pae = np.asarray(payload["pae"], dtype=np.float32)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"PAE in {confidences_json} is not square (shape {pae.shape})")

    res_ids = payload.get("token_res_ids")
    chain_ids = payload.get("token_chain_ids")
    if res_ids is None:
        residue_ids = np.arange(1, pae.shape[0] + 1, dtype=np.int32)
    else:
        residue_ids = np.asarray(res_ids, dtype=np.int32)
    if chain_ids is not None and len(set(chain_ids)) > 1:
        print(
            f"  [warn] {job_id}: multi-chain ({sorted(set(chain_ids))}); token_res_ids "
            "may repeat across chains. interdomain_pae assumes single-chain TPS designs."
        )
    if residue_ids.shape[0] != pae.shape[0]:
        raise ValueError(
            f"token_res_ids length {residue_ids.shape[0]} != PAE dim {pae.shape[0]} "
            f"in {confidences_json}"
        )

    ptm, iptm = _load_ptm_iptm(confidences_json, job_id)

    os.makedirs(os.path.dirname(os.path.abspath(out_npz)), exist_ok=True)
    np.savez_compressed(
        out_npz,
        pae=pae,
        residue_ids=residue_ids,
        n_residues=np.int64(pae.shape[0]),
        source="alphafold3",
        ptm=np.float32(ptm),
        iptm=np.float32(iptm),
    )
    return int(pae.shape[0])


def extract_af_output(af_output: Path, pae_dir: str, *, skip_existing: bool = True) -> int:
    """Write a ``<job>_pae.npz`` for every AF3 job under `af_output`. Returns count."""
    os.makedirs(pae_dir, exist_ok=True)
    n = 0
    for entry in sorted(os.listdir(af_output)):
        conf = af_output / entry / f"{entry}_confidences.json"
        if not conf.is_file():
            continue
        out_npz = os.path.join(pae_dir, f"{entry}_pae.npz")
        if skip_existing and os.path.isfile(out_npz):
            print(f"[skip] {entry}: {out_npz} already exists")
            n += 1
            continue
        try:
            L = confidences_to_npz(str(conf), out_npz, job_id=entry)
        except Exception as exc:  # noqa: BLE001 - keep extracting the rest
            print(f"  [warn] {entry}: {exc}")
            continue
        print(f"[ok] {entry}: PAE ({L}x{L}) -> {out_npz}")
        n += 1
    print(f"Wrote/kept {n} PAE npz file(s) in {pae_dir}")
    return n


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract AlphaFold3 PAE matrices from an af_output tree into the "
        "shared <ID>_pae.npz schema (pae, residue_ids, n_residues, source) consumed "
        "by interdomain_pae.py. ID = af_output job subfolder name."
    )
    parser.add_argument("--af_output", default=None, help="AF3 af_output directory.")
    parser.add_argument(
        "--structs_dir",
        default=None,
        help="Alternatively, a structs dir whose sibling/child af_output/ is used.",
    )
    parser.add_argument(
        "--pae_dir", required=True, help="Output directory for <ID>_pae.npz files."
    )
    parser.add_argument(
        "--no-skip_existing",
        dest="skip_existing",
        action="store_false",
        help="Re-extract even if <ID>_pae.npz exists (default: skip existing).",
    )
    args = parser.parse_args()
    af_output = _resolve_af_output(args.af_output, args.structs_dir)
    print(f"AF3 af_output: {af_output}")
    extract_af_output(af_output, args.pae_dir, skip_existing=args.skip_existing)


if __name__ == "__main__":
    main()
