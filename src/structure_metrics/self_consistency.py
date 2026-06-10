from __future__ import annotations

# Self-consistency scRMSD — the headline designability metric.
#
# For each design backbone:
#   1. Sample N sequences from the backbone with ProteinMPNN (default N=8).
#   2. Refold each sampled sequence with ESMFold (reuses src/esmfold).
#   3. Cα-align each refold back to the ORIGINAL design structure (Biopython
#      Superimposer) and record the RMSD.
# Report sc_rmsd_min (best of N) and sc_rmsd_mean. A design is "self-consistent"
# / designable when sc_rmsd_min < ~2 Angstrom: there exists a sequence that both
# ProteinMPNN likes for this fold AND that ESMFold folds back to the same shape.
#
# This is GPU + slow (N folds per structure). N (`--num_seqs`) and the structure
# subset (`--ids` / `--limit`) are configurable so it can be validated on 1-2
# structures cheaply. Runs in the `esmfold` conda env (has torch + transformers +
# Biopython); ProteinMPNN is shelled out with the SAME python so no extra env.
#
# Output (CSV keyed by ID): sc_rmsd_min, sc_rmsd_mean, n_samples (folds that
# succeeded). ID = structure filename stem (matches plddt / proteinmpnn_score).

import glob
import os
import subprocess
import sys
import tempfile
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select, Superimposer
from Bio.PDB.PDBExceptions import PDBConstructionWarning

SRC_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SRC_DIR.parent
PROTEINMPNN_DIR = REPO_ROOT / "vendor" / "ProteinMPNN"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

COLUMNS = ["ID", "sc_rmsd_min", "sc_rmsd_mean", "n_samples"]

_PDB_PARSER = PDBParser(QUIET=True)
_CIF_PARSER = MMCIFParser(QUIET=True)


def _parser_for(path: str):
    return _CIF_PARSER if path.lower().endswith((".cif", ".mmcif")) else _PDB_PARSER


def _chain_ids(structure_path: str) -> List[str]:
    """Chain IDs (first model) that contain at least one standard polymer residue."""
    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    model = next(iter(structure))
    ids = []
    for chain in model:
        if any(res.id[0] == " " and "CA" in res for res in chain):
            ids.append(chain.id)
    return ids


class _ChainSelect(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):  # noqa: N802 (Biopython API)
        return 1 if chain.id == self.chain_id else 0

    def accept_residue(self, residue):  # noqa: N802
        return 1 if residue.id[0] == " " else 0


def _write_single_chain(structure_path: str, chain_id: str, out_path: str) -> None:
    """Write only `chain_id` (standard residues) of `structure_path` as a PDB.
    scRMSD is a single-chain (monomer) designability metric; multimer inputs are
    reduced to one design chain so ProteinMPNN designs / ESMFold folds / RMSD
    compares the SAME single chain (otherwise a complex vs a folded-as-monomer
    concatenation gives meaningless RMSD)."""
    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, _ChainSelect(chain_id))


def _ca_atoms(structure_path: str):
    """Ordered list of Cα atoms (first model, standard residues, all chains)."""
    parser = _parser_for(structure_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        structure = parser.get_structure("s", structure_path)
    model = next(iter(structure))
    cas = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if "CA" in residue:
                cas.append(residue["CA"])
    return cas


def _ca_rmsd(ref_path: str, mobile_path: str) -> float:
    """Cα-RMSD of `mobile` superimposed onto `ref`. Requires equal residue counts
    (true for a refold of a backbone-derived sequence); if they differ we align
    the leading min(len) atoms and warn."""
    ref = _ca_atoms(ref_path)
    mob = _ca_atoms(mobile_path)
    if not ref or not mob:
        return float("nan")
    if len(ref) != len(mob):
        n = min(len(ref), len(mob))
        print(f"  [warn] residue count mismatch ({len(ref)} vs {len(mob)}); aligning first {n}")
        ref, mob = ref[:n], mob[:n]
    sup = Superimposer()
    sup.set_atoms(ref, mob)
    return float(sup.rms)


def _collect_structures(structs_dir: str):
    """ID -> structure file. Mirrors plddt.py's af3-vs-flat detection."""
    af3: Dict[str, str] = {}
    try:
        entries = sorted(os.listdir(structs_dir))
    except FileNotFoundError:
        entries = []
    for entry in entries:
        sub = os.path.join(structs_dir, entry)
        model = os.path.join(sub, entry + "_model.cif")
        if os.path.isdir(sub) and os.path.isfile(model):
            af3[entry] = model
    if af3:
        return OrderedDict(sorted(af3.items())), "af3"

    chosen: Dict[str, str] = {}
    for ext in (".mmcif", ".cif", ".pdb"):
        for path in sorted(glob.glob(os.path.join(structs_dir, f"*{ext}"))):
            stem = os.path.splitext(os.path.basename(path))[0]
            chosen[stem] = path
    return OrderedDict(sorted(chosen.items())), "flat"


def _default_save_path(structs_dir: str) -> str:
    d = structs_dir.rstrip(os.sep)
    return os.path.join(os.path.dirname(d), os.path.basename(d) + "_self_consistency.csv")


def _sample_sequences(
    pdb_path: str,
    out_folder: str,
    *,
    num_seqs: int,
    sampling_temp: float,
    model_name: str,
    seed: int,
    pdb_path_chains: Optional[str] = None,
    python_exe: str = sys.executable,
) -> List[str]:
    """Run ProteinMPNN in sampling mode; return the N sampled sequences (chains
    joined, '/' chain separators dropped -> single-chain folding input)."""
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    cmd = [
        python_exe,
        str(PROTEINMPNN_DIR / "protein_mpnn_run.py"),
        "--pdb_path", pdb_path,
        "--out_folder", out_folder,
        "--num_seq_per_target", str(num_seqs),
        "--sampling_temp", str(sampling_temp),
        "--model_name", model_name,
        "--seed", str(seed),
        "--batch_size", "1",
    ]
    if pdb_path_chains:
        cmd += ["--pdb_path_chains", pdb_path_chains]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"ProteinMPNN sampling failed on {pdb_path}:\n{proc.stdout}\n{proc.stderr}"
        )
    fa_path = os.path.join(out_folder, "seqs", stem + ".fa")
    if not os.path.isfile(fa_path):
        raise RuntimeError(f"ProteinMPNN produced no seqs file {fa_path}\n{proc.stdout}")

    # The .fa has the native sequence FIRST (T=...), then N sampled sequences.
    seqs: List[str] = []
    with open(fa_path) as fh:
        header = None
        cur = []
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append(("".join(cur), header))
                header = line
                cur = []
            else:
                cur.append(line)
        if header is not None:
            seqs.append(("".join(cur), header))
    # Drop the native (first record); keep sampled ones; join multi-chain on '/'.
    sampled = [s.replace("/", "") for s, h in seqs[1:]]
    return sampled


def self_consistency_for_structure(
    pdb_path: str,
    *,
    fold_fn,
    num_seqs: int,
    sampling_temp: float,
    model_name: str,
    seed: int,
    workdir: str,
    chain: Optional[str] = None,
) -> Dict[str, float]:
    """scRMSD for a single design: sample -> refold -> Cα-RMSD to original.

    scRMSD is a single-chain (monomer) designability metric. If the structure has
    multiple chains we restrict to ONE design chain (`chain`, default = first):
    ProteinMPNN designs only that chain, ESMFold folds the single sampled chain,
    and the Cα-RMSD compares to that same chain. (Comparing a multi-chain complex
    to a folded-as-monomer concatenation otherwise gives a meaningless RMSD.)"""
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    mpnn_out = os.path.join(workdir, "mpnn", stem)
    os.makedirs(mpnn_out, exist_ok=True)

    chains = _chain_ids(pdb_path)
    if not chains:
        return {"sc_rmsd_min": float("nan"), "sc_rmsd_mean": float("nan"), "n_samples": 0}
    design_chain = chain if (chain and chain in chains) else chains[0]
    if len(chains) > 1:
        print(f"  [note] {stem}: {len(chains)} chains {chains}; scoring single chain "
              f"'{design_chain}' (monomer designability).")
        ref_path = os.path.join(workdir, f"{stem}_chain_{design_chain}.pdb")
        _write_single_chain(pdb_path, design_chain, ref_path)
        pdb_chains_arg = design_chain
    else:
        ref_path = pdb_path
        pdb_chains_arg = None

    sequences = _sample_sequences(
        ref_path, mpnn_out, num_seqs=num_seqs, sampling_temp=sampling_temp,
        model_name=model_name, seed=seed, pdb_path_chains=pdb_chains_arg,
    )
    rmsds: List[float] = []
    for k, seq in enumerate(sequences):
        refold_path = os.path.join(workdir, f"{stem}_refold_{k}.pdb")
        try:
            pdb_str = fold_fn(seq)
            with open(refold_path, "w") as fh:
                fh.write(pdb_str)
            rmsd = _ca_rmsd(ref_path, refold_path)
        except Exception as exc:  # noqa: BLE001
            print(f"  [warn] refold/RMSD failed for {stem} sample {k}: {exc}")
            continue
        rmsds.append(rmsd)
        print(f"  {stem} sample {k}: scRMSD = {rmsd:.3f} A")
    if not rmsds:
        return {"sc_rmsd_min": float("nan"), "sc_rmsd_mean": float("nan"), "n_samples": 0}
    arr = np.asarray(rmsds, dtype=float)
    return {
        "sc_rmsd_min": float(np.nanmin(arr)),
        "sc_rmsd_mean": float(np.nanmean(arr)),
        "n_samples": int(np.isfinite(arr).sum()),
    }


def self_consistency_dir(
    structs_dir: str,
    *,
    save_path: Optional[str] = None,
    num_seqs: int = 8,
    sampling_temp: float = 0.1,
    model_name: str = "v_48_020",
    seed: int = 0,
    ids: Optional[List[str]] = None,
    limit: Optional[int] = None,
    device: Optional[str] = None,
    chain: Optional[str] = None,
) -> pd.DataFrame:
    """Self-consistency scRMSD for structures in a dir. Loads ESMFold ONCE and
    reuses it across all structures. `ids`/`limit` restrict the structure subset
    (validate on 1-2 first — this is GPU + slow)."""
    structures, mode = _collect_structures(structs_dir)
    if not structures:
        raise ValueError(
            f"No structures found in {structs_dir} (expected AF3 af_output or flat .pdb/.cif)."
        )
    if ids:
        wanted = set(ids)
        structures = OrderedDict((k, v) for k, v in structures.items() if k in wanted)
        missing = wanted - set(structures)
        if missing:
            print(f"  [warn] requested IDs not found: {sorted(missing)}")
    if limit is not None:
        structures = OrderedDict(list(structures.items())[:limit])
    if not structures:
        raise ValueError("No structures left to score after --ids/--limit filtering.")
    print(f"Detected {mode} layout: scoring {len(structures)} structure(s) in {structs_dir}")

    # Load ESMFold once (reuse src/esmfold).
    from esmfold.esmfold import _load_model, fold_sequence

    tokenizer, model, device = _load_model(device)

    def fold_fn(seq: str) -> str:
        return fold_sequence(model, tokenizer, device, seq)

    rows: List[Dict[str, float]] = []
    n = len(structures)
    with tempfile.TemporaryDirectory(prefix="self_consistency_") as workdir:
        for i, (stem, path) in enumerate(structures.items(), start=1):
            print(f"[{i}/{n}] {stem}: sampling {num_seqs} seqs + refolding ...")
            try:
                stats = self_consistency_for_structure(
                    path, fold_fn=fold_fn, num_seqs=num_seqs, sampling_temp=sampling_temp,
                    model_name=model_name, seed=seed, workdir=workdir, chain=chain,
                )
            except Exception as exc:  # noqa: BLE001
                print(f"  [warn] failed on {stem}: {exc}")
                stats = {"sc_rmsd_min": float("nan"), "sc_rmsd_mean": float("nan"), "n_samples": 0}
            stats["ID"] = stem
            rows.append(stats)
            print(f"  -> {stem}: sc_rmsd_min={stats['sc_rmsd_min']}, "
                  f"sc_rmsd_mean={stats['sc_rmsd_mean']}, n={stats['n_samples']}")

    df = pd.DataFrame(rows)[COLUMNS].sort_values("ID").reset_index(drop=True)

    if save_path is None:
        save_path = _default_save_path(structs_dir)
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
