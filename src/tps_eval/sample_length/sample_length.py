import argparse

from Bio import SeqIO
import numpy as np
import csv

def get_sampler(fasta_path:str, return_counts:bool=False):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    lengths = [len(record.seq) for record in records]
    unique_lengths, counts = np.unique(lengths, return_counts=True)

    if return_counts:
        def sample_lengths(n):
            """
            Example usage:
                sampled_lengths, counts = sample_lengths(n)
            """
            sampled_lengths = np.random.choice(unique_lengths, size=n, p=counts/counts.sum())
            unique, sampled_counts = np.unique(sampled_lengths, return_counts=True)
            return list(zip(unique, sampled_counts))
    else:
        def sample_lengths(n):
            """
            Example usage:
                sampled_lengths = sample_lengths(n)
            """
            return np.random.choice(unique_lengths, size=n, p=counts/counts.sum())

    return sample_lengths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_path", type=str)
    parser.add_argument("--num_seqs", type=int)
    parser.add_argument("--save_to", type=str)
    parser.add_argument("--return_counts", default=False, action=argparse.BooleanOptionalAction)
    args = parser.parse_args()

    sampler = get_sampler(args.fasta_path, return_counts=args.return_counts)
    
    if args.return_counts:
        unique_lengths = sampler(args.num_seqs)
    else:
        lengths = sampler(args.num_seqs)
    
    if args.return_counts:
        with open(args.save_to, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['length', 'count'])
            for length, count in unique_lengths:
                writer.writerow([length, count])
    else:   
        with open(args.save_to, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['length'])
            for length in lengths:
                writer.writerow([length])


if __name__ == "__main__":
    main()
