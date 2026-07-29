import os
import argparse
from pathlib import Path
import logging
import subprocess
from uuid import uuid4
from shutil import rmtree

import pandas as pd  # type: ignore

from tps_eval.pandas_compat import group_idxmax_skipna

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
    
    foldseek_command = f'foldseek easy-search {args.detected_domain_structures_root} {args.known_domain_structures_root} {tsv_path} {tmp_path} --max-seqs 5000 -e 1 -s 10 --exhaustive-search -v 3 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,lddt'.split()
    process = subprocess.Popen(foldseek_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in process.stdout:
        print(line, end="")
    process.stdout.close()
    return_code = process.wait()

    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, foldseek_command)
    print("Foldseek finished successfully.")
    
    # Create final output CSV
    df_foldseek = pd.read_csv(tsv_path, sep='\t', header=None,
                              names=['query', 'target', 'fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend',
                                     'tstart', 'tend', 'evalue', 'bits', 'alntmscore', 'qtmscore', 'ttmscore', 'lddt'])
    # Foldseek's score columns can be read back as object dtype (e.g. an lddt cell
    # with a stray non-numeric token), which makes the groupby idxmax below raise
    # "'>' not supported between instances of 'str' and 'float'". Coerce the numeric
    # columns so idxmax works (non-numeric -> NaN, skipped by idxmax).
    for _num_col in ['fident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend',
                     'tstart', 'tend', 'evalue', 'bits', 'alntmscore', 'qtmscore',
                     'ttmscore', 'lddt']:
        df_foldseek[_num_col] = pd.to_numeric(df_foldseek[_num_col], errors='coerce')
    df_foldseek.to_csv(tsv_path.with_suffix('.csv'), index=False)
    dfg = df_foldseek.groupby('query')

    def _val_at(idx, col):
        # A query whose entire score column is NaN yields a NaN idxmax (see
        # group_idxmax_skipna -- pandas 3 would otherwise RAISE on such a group); guard
        # the lookup so it produces a NaN cell instead of raising KeyError on .loc[NaN].
        return df_foldseek.loc[idx, col] if pd.notna(idx) else float("nan")

    best_scores = pd.DataFrame({
        'query': dfg.groups.keys(),
        'max_alntmscore': group_idxmax_skipna(dfg, 'alntmscore').map(lambda idx: _val_at(idx, 'alntmscore')),
        'max_alntmscore_target': group_idxmax_skipna(dfg, 'alntmscore').map(lambda idx: _val_at(idx, 'target')),
        'max_qtmscore': group_idxmax_skipna(dfg, 'qtmscore').map(lambda idx: _val_at(idx, 'qtmscore')),
        'max_qtmscore_target': group_idxmax_skipna(dfg, 'qtmscore').map(lambda idx: _val_at(idx, 'target')),
        'max_ttmscore': group_idxmax_skipna(dfg, 'ttmscore').map(lambda idx: _val_at(idx, 'ttmscore')),
        'max_ttmscore_target': group_idxmax_skipna(dfg, 'ttmscore').map(lambda idx: _val_at(idx, 'target')),
        'max_lddt': group_idxmax_skipna(dfg, 'lddt').map(lambda idx: _val_at(idx, 'lddt')),
        'max_lddt_target': group_idxmax_skipna(dfg, 'lddt').map(lambda idx: _val_at(idx, 'target')),
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
