from __future__ import annotations

import argparse

from self_consistency import self_consistency_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Self-consistency scRMSD (designability) for design structures. "
        "For each backbone: sample N sequences with ProteinMPNN, refold each with "
        "ESMFold, and report the min/mean Cα-RMSD of a refold back to the original. "
        "Designable ~ sc_rmsd_min < 2 A. GPU + slow; restrict with --limit/--ids."
    )
    parser.add_argument(
        "structs_dir",
        help="AF3 af_output dir (<job>/<job>_model.cif) OR flat dir of .pdb/.cif; "
        "ID = job name / filename stem.",
    )
    parser.add_argument(
        "--save_path", default=None,
        help="Output CSV path (default: <structs_dir>_self_consistency.csv).",
    )
    parser.add_argument("--num_seqs", type=int, default=8,
                        help="ProteinMPNN sequences sampled per backbone (default 8).")
    parser.add_argument("--sampling_temp", type=float, default=0.1,
                        help="ProteinMPNN sampling temperature (default 0.1).")
    parser.add_argument("--model_name", default="v_48_020",
                        help="ProteinMPNN model name (default v_48_020).")
    parser.add_argument("--seed", type=int, default=0, help="Random seed (0 = random).")
    parser.add_argument("--ids", nargs="+", default=None,
                        help="Restrict to these structure IDs (validate on 1-2 first).")
    parser.add_argument("--limit", type=int, default=None,
                        help="Score only the first N structures (cheap validation).")
    parser.add_argument("--device", default=None, help="Torch device for ESMFold (cuda/cpu).")
    parser.add_argument("--chain", default=None,
                        help="For multi-chain structures, the single design chain to score "
                        "(default: first chain). scRMSD is a single-chain designability metric.")
    args = parser.parse_args()

    self_consistency_dir(
        args.structs_dir,
        save_path=args.save_path,
        num_seqs=args.num_seqs,
        sampling_temp=args.sampling_temp,
        model_name=args.model_name,
        seed=args.seed,
        ids=args.ids,
        limit=args.limit,
        device=args.device,
        chain=args.chain,
    )


if __name__ == "__main__":
    main()
