from __future__ import annotations

import argparse

from esmfold import fold_fasta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Predict structures with ESMFold (facebook/esmfold_v1) for "
        "every sequence in a FASTA file, writing one <ID>.pdb per record into a "
        "structs dir. The PDB B-factor field holds per-residue pLDDT (0-100, like "
        "AlphaFold), so the output dir is consumed unchanged by run_plddt.sh and "
        "run_structural_identity.sh (ID = FASTA record id = filename stem)."
    )
    parser.add_argument(
        "fasta_path",
        help="Path to the FASTA file of sequences to fold.",
    )
    parser.add_argument(
        "--save_dir",
        required=True,
        help="Directory to write <ID>.pdb structures into (created if missing). "
        "Point run_plddt.sh --structs_dir at this directory afterwards.",
    )
    parser.add_argument(
        "--no-skip_existing",
        dest="skip_existing",
        action="store_false",
        help="Re-fold sequences even if <ID>.pdb already exists (default: skip existing).",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=None,
        help="Force ESMFold trunk chunk size (lower = less GPU memory, slower). "
        "Default: auto (no chunking for short seqs; chunked for long seqs).",
    )
    parser.add_argument(
        "--device",
        default=None,
        help="Torch device ('cuda'/'cpu'). Default: cuda if available, else cpu.",
    )
    args = parser.parse_args()
    fold_fasta(
        args.fasta_path,
        args.save_dir,
        skip_existing=args.skip_existing,
        chunk_size=args.chunk_size,
        device=args.device,
    )


if __name__ == "__main__":
    main()
