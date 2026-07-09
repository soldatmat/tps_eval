from __future__ import annotations

import argparse

from tps_eval.sequence_metrics.tmprot import score_fasta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="TmProt melting-temperature (Tm) prediction. Writes "
        "<input>_tmprot.csv keyed by ID with the RAW predicted Tm (degC) in the "
        "`tm` column (one row per input FASTA record; NaN where TmProt skipped a "
        "sequence it cannot score)."
    )
    parser.add_argument("--fasta_path", required=True, help="Path to the FASTA file.")
    parser.add_argument(
        "--device",
        default=None,
        help="Optional device override: 'cpu' hides the GPU; 'cuda'/omitted "
        "leaves TmProt's automatic cuda-if-available selection.",
    )
    parser.add_argument(
        "--out_suffix",
        default="tmprot",
        help="Output filename suffix (default: tmprot -> <input>_tmprot.csv).",
    )
    args = parser.parse_args()

    score_fasta(
        args.fasta_path,
        out_suffix=args.out_suffix,
        device=args.device,
    )


if __name__ == "__main__":
    main()
