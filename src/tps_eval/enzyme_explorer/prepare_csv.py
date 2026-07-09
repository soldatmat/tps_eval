import pandas as pd
from Bio import SeqIO
import argparse
import os

def fasta_to_csv(fasta_path, csv_path):
    records = [{'ID': rec.id, 'sequence': str(rec.seq)} for rec in SeqIO.parse(fasta_path, "fasta")]
    df = pd.DataFrame(records)
    df.to_csv(csv_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert FASTA to CSV.")
    parser.add_argument("--fasta_path", required=True, help="Path to the input FASTA file.")
    parser.add_argument("--csv_path", required=False, help="Path to the output CSV file.")
    args = parser.parse_args()

    fasta_path = args.fasta_path
    if args.csv_path:
        csv_path = args.csv_path
    else:
        csv_path = os.path.splitext(fasta_path)[0] + ".csv"

    fasta_to_csv(fasta_path, csv_path)
