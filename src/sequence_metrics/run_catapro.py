from __future__ import annotations

import argparse

from catapro import DEFAULT_OUT_SUFFIX, score_fasta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CataPro enzyme-kinetics prediction (kcat, Km, kcat/Km) for every "
        "sequence in a FASTA file, against one substrate. Writes a CSV keyed by ID "
        "(catapro_kcat [s^-1], catapro_km [mM], catapro_kcat_km [s^-1 mM^-1], "
        "catapro_substrate). The substrate SMILES is resolved from --target_substrate "
        "via alphafold.cofold_substrates.SUBSTRATE_SMILES (or given directly with "
        "--smiles); a substrate with no known SMILES yields NaN rows."
    )
    parser.add_argument("fasta_path", help="FASTA file of sequences to score.")
    parser.add_argument(
        "--target_substrate",
        default=None,
        help="Substrate code (e.g. GPP, FPP, GGPP, GFPP) resolved to a SMILES.",
    )
    parser.add_argument(
        "--smiles",
        default=None,
        help="Explicit substrate SMILES (overrides --target_substrate lookup).",
    )
    parser.add_argument(
        "--model_dpath",
        default=None,
        help="CataPro model directory (default: vendor/CataPro/models).",
    )
    parser.add_argument(
        "--device",
        default="cuda",
        help="Torch device for CataPro/ProtT5/MolT5 (default: cuda; use cpu on CPU nodes).",
    )
    parser.add_argument("--batch_size", type=int, default=64, help="CataPro batch size.")
    parser.add_argument(
        "--out_suffix",
        default=DEFAULT_OUT_SUFFIX,
        help=f"Output filename suffix -> <fasta_stem>_<suffix>.csv (default: {DEFAULT_OUT_SUFFIX}; "
        "the MARTS band panel uses e.g. catapro_FPP).",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <fasta_stem>_<out_suffix>.csv beside the FASTA).",
    )
    args = parser.parse_args()

    score_fasta(
        args.fasta_path,
        args.target_substrate,
        smiles=args.smiles,
        model_dpath=args.model_dpath,
        device=args.device,
        batch_size=args.batch_size,
        out_suffix=args.out_suffix,
        save_path=args.save_path,
    )


if __name__ == "__main__":
    main()
