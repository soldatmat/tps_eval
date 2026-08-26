from __future__ import annotations

import argparse
import json

from tps_eval.enzyme_explorer.enzyme_explorer_structure import (
    default_save_path,
    reshape_to_csv,
    summarize,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reshape EnzymeExplorer's `predict_with_structures` output (two CSVs: "
        "predictions_plm_domains.csv + predictions_plm_only_fallback.csv) into ONE CSV "
        "keyed by ID with one row per design. A design whose structure yielded NO TPS "
        "structural domain is absent from the domains CSV; it gets ee_struct_status="
        "'no_domains' and a NaN isTPS_struct_* score, so any >=threshold test fails it. "
        "The PLM-only fallback score is carried under isTPS_fallback_* for DIAGNOSTICS "
        "ONLY -- it is a different classifier and must never substitute for the structure "
        "score. Pure pandas; no EnzymeExplorer import needed."
    )
    parser.add_argument(
        "output_dir",
        help="The --output-dir that `predict_with_structures` wrote its two CSVs into.",
    )
    parser.add_argument(
        "--structs_dir",
        required=True,
        help="Structures directory the prediction ran on; sets the default output path "
        "(<structs_dir>_enzyme_explorer_structure.csv) and is the fallback ID universe.",
    )
    parser.add_argument(
        "--fasta_path",
        default=None,
        help="Input FASTA -- the authoritative design universe (preferred: it includes "
        "designs that never folded, which SHOULD appear as failures).",
    )
    parser.add_argument(
        "--sequences_csv_path",
        default=None,
        help="Input sequences CSV, as an alternative design universe to --fasta_path.",
    )
    parser.add_argument(
        "--csv_id_column",
        default="id",
        help="ID column in --sequences_csv_path (default: id).",
    )
    parser.add_argument(
        "--save_path",
        default=None,
        help="Output CSV path (default: <structs_dir>_enzyme_explorer_structure.csv).",
    )
    parser.add_argument(
        "--id_map_csv",
        default=None,
        help="CSV with columns sanitized_id,original_id -- restores the original design "
        "IDs after a run that had to sanitize them. EnzymeExplorer rejects any structure "
        "filename stem containing a character outside [a-zA-Z0-9_] (it raises on the first "
        "offender and aborts the whole directory), so arms whose IDs contain dots or "
        "hyphens run through a sanitized symlink farm.",
    )
    args = parser.parse_args()

    path = reshape_to_csv(
        args.output_dir,
        structs_dir=args.structs_dir,
        fasta_path=args.fasta_path,
        sequences_csv_path=args.sequences_csv_path,
        csv_id_column=args.csv_id_column,
        save_path=args.save_path or default_save_path(args.structs_dir),
        id_map_csv=args.id_map_csv,
    )
    import pandas as pd

    print(f"Wrote {path}")
    print(json.dumps(summarize(pd.read_csv(path)), indent=2))


if __name__ == "__main__":
    main()
