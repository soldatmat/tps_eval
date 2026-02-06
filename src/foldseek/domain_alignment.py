import os
import argparse
from pathlib import Path
import logging
import subprocess
from uuid import uuid4
from shutil import rmtree

import pandas as pd  # type: ignore

logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def parse_args() -> argparse.Namespace:
    """
    This function parses arguments
    :return: current argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="A script to compare detected TPS domains to the known ones")
    parser.add_argument("--known_domain_structures_root", help="A directory containing structures of known domains", type=str)
    parser.add_argument("--detected_domain_structures_root",help="A path to new detected domain structures",type=str)
    parser.add_argument("--output_root", type=str, required=True, help="Path to output CSV file.")
    parser.add_argument("--store_intermediate_results", action="store_true", help="Flag to keep files with intermediate results.", default=False)
    parser.add_argument("--random_run_id", action="store_true", default=False, help="Flag to add random uuid4 to output files to avoid overwriting.")
    return parser.parse_args()


def main(args: argparse.Namespace):
    run_id = None
    if args.random_run_id:
        run_id = str(uuid4())

    # Run Foldseek
    output_root = Path(args.output_root)
    if not output_root.exists():
        output_root.mkdir()
    tsv_path = output_root / f'domain_alignments{f"_{run_id}" if args.random_run_id else ""}.tsv'
    tmp_path = output_root / f'tmp{f"_{run_id}" if args.random_run_id else ""}'
    foldseek_comparison_output = subprocess.check_output(
        f'foldseek easy-search {args.detected_domain_structures_root} {args.known_domain_structures_root} {tsv_path} {tmp_path} --max-seqs 5000 -e 1 -s 10 --exhaustive-search --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,lddt'.split())
    logger.info(f'Foldseek output: {foldseek_comparison_output.decode("utf-8")}')
    
    # Create final output CSV
    df_foldseek = pd.read_csv(tsv_path, sep='\t', header=None,
                              names=['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend',
                                     'tstart', 'tend', 'evalue', 'bits', 'alntmscore', 'qtmscore', 'ttmscore', 'lddt'])
    df_foldseek.to_csv(tsv_path.with_suffix('.csv'), index=False)
    dfg = df_foldseek.groupby('query')
    dfg.groups.keys(), dfg['alntmscore'].max(), dfg['alntmscore'].idxmax()
    best_scores = pd.DataFrame({
        'query': dfg.groups.keys(),

        # 'max_alntmscore_idx': dfg['alntmscore'].idxmax().values,
        'max_alntmscore': dfg['alntmscore'].idxmax().map(lambda idx: df_foldseek.loc[idx, 'alntmscore']),
        'max_alntmscore_target': dfg['alntmscore'].idxmax().map(lambda idx: df_foldseek.loc[idx, 'target']),

        # 'max_qtmscore_idx': dfg['qtmscore'].idxmax().values,
        'max_qtmscore': dfg['qtmscore'].idxmax().map(lambda idx: df_foldseek.loc[idx, 'qtmscore']),
        'max_qtmscore_target': dfg['qtmscore'].idxmax().map(lambda idx: df_foldseek.loc[idx, 'target']),

        # 'max_ttmscore_idx': dfg['ttmscore'].idxmax().values,
        'max_ttmscore': dfg['ttmscore'].idxmax().map(lambda idx: df_foldseek.loc[idx, 'ttmscore']),
        'max_ttmscore_target': dfg['ttmscore'].idxmax().map(lambda idx: df_foldseek.loc[idx, 'target']),

        # 'max_lddt_idx': dfg['lddt'].idxmax().values,
        'max_lddt': dfg['lddt'].idxmax().map(lambda idx: df_foldseek.loc[idx, 'lddt']),
        'max_lddt_target': dfg['lddt'].idxmax().map(lambda idx: df_foldseek.loc[idx, 'target']),
    })
    csv_path = output_root / "domain_alignment_scores.csv"
    if args.random_run_id:
        csv_path = output_root / (csv_path.stem + f'_{run_id}' + csv_path.suffix)
    best_scores.to_csv(csv_path, index=False)

    # Clean up intermediate results
    if not args.store_intermediate_results:
        os.remove(tsv_path)
        os.remove(tsv_path.with_suffix('.csv'))
    rmtree(tmp_path)


if __name__ == "__main__":
    args = parse_args()
    main(args)
