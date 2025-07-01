from Bio import SeqIO
import numpy as np

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
