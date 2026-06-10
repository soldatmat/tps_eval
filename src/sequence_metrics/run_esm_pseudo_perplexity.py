from __future__ import annotations

import argparse

from esm_pseudo_perplexity import DEFAULT_MODEL, compute_pseudo_perplexity


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ESM pseudo-perplexity (naturalness) for every sequence in a "
        "FASTA file. Lower = more natural / in-distribution under ESM's masked LM. "
        "Writes a CSV keyed by ID (esm_pseudo_perplexity, esm_mean_pll). Uses the "
        "same ESM-1b model as the embedding tool."
    )
    parser.add_argument("fasta_path", help="FASTA file of sequences to score.")
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <fasta>_esm_pseudo_perplexity.csv).",
    )
    parser.add_argument(
        "--model_location",
        default=DEFAULT_MODEL,
        help=f"ESM model name or path (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--method",
        choices=["swoop", "masked"],
        default="swoop",
        help="'swoop' = One Fell Swoop single-pass approximation (fast, default); "
        "'masked' = exact masked-marginal pseudo-perplexity (O(L) forwards/seq, slow).",
    )
    parser.add_argument("--toks_per_batch", type=int, default=4096, help="Max batch token budget (swoop).")
    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available.")
    args = parser.parse_args()

    compute_pseudo_perplexity(
        args.fasta_path,
        save_path=args.save_path,
        model_location=args.model_location,
        method=args.method,
        toks_per_batch=args.toks_per_batch,
        nogpu=args.nogpu,
    )


if __name__ == "__main__":
    main()
