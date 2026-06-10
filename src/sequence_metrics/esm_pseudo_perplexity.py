from __future__ import annotations

# ESM pseudo-perplexity (naturalness) metric.
#
# For each sequence we score how "in-distribution" / natural it is under ESM's
# masked language model. Lower pseudo-perplexity == the sequence's residues are,
# on average, well-predicted by ESM from their context == more natural.
#
# Two estimators of the per-residue pseudo-log-likelihood (PLL) are available:
#
#  * "swoop" (default, fast) — the "One Fell Swoop" single-forward-pass
#    approximation (Salazar et al. 2020 PLL approximated by NOT masking; cf.
#    "Single Sequence Prediction over Reasoning ... ESM" usage): run the model
#    ONCE on the unmasked sequence and read off log p(x_i | x_{!=i}) ~= the
#    logit at the true residue at each position. This skips the O(L) per-residue
#    masking loop, so a whole FASTA scores in one batched pass. It slightly
#    underestimates perplexity (the model peeks at the residue it is scoring) but
#    rank-correlates strongly with the exact masked-marginal score and is the
#    standard cheap naturalness proxy.
#
#  * "masked" (exact, slow) — true masked marginals: for each position, mask it,
#    forward, read log p(true | context). O(L) forwards per sequence. Use for a
#    handful of sequences when you need the exact pseudo-perplexity.
#
# Output (CSV keyed by ID): esm_mean_pll (mean per-residue log-likelihood, <=0;
# higher = more natural) and esm_pseudo_perplexity = exp(-esm_mean_pll)
# (>=1; LOWER = more natural). Same model as src/esm/extract_embeddings.py
# (ESM-1b, esm1b_t33_650M_UR50S) so the naturalness score is consistent with the
# embedding-based metrics.

import os
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import torch

from esm import pretrained

SRC_DIR = Path(__file__).resolve().parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from data.sequences import load_fasta_sequences, separate_identifiers  # noqa: E402

# Same model the embedding tool loads, so naturalness is consistent across metrics.
DEFAULT_MODEL = "esm1b_t33_650M_UR50S"
# ESM-1b was trained on sequences truncated to 1022 residues (1024 with BOS/EOS).
TRUNCATION_SEQ_LENGTH = 1022

COLUMNS = ["ID", "esm_pseudo_perplexity", "esm_mean_pll", "n_residues"]


def _per_position_log_probs(logits: torch.Tensor) -> torch.Tensor:
    """Log-softmax over the vocabulary axis (last). logits: (..., vocab)."""
    return torch.log_softmax(logits, dim=-1)


def _score_swoop(model, alphabet, device, sequences: List[str], toks_per_batch: int) -> List[dict]:
    """One Fell Swoop: a single (batched) unmasked forward pass per sequence,
    reading log p(x_i | full context) at each residue's true token."""
    batch_converter = alphabet.get_batch_converter(TRUNCATION_SEQ_LENGTH)
    results: List[dict] = []

    # Batch by token budget (mirrors extract_embeddings' FastaBatchedDataset idea
    # without needing labels in the batch).
    idx_order = sorted(range(len(sequences)), key=lambda i: len(sequences[i]))
    batches: List[List[int]] = []
    cur: List[int] = []
    cur_max = 0
    for i in idx_order:
        L = min(len(sequences[i]), TRUNCATION_SEQ_LENGTH) + 2  # + BOS/EOS
        new_max = max(cur_max, L)
        if cur and new_max * (len(cur) + 1) > toks_per_batch:
            batches.append(cur)
            cur, cur_max = [], 0
            new_max = L
        cur.append(i)
        cur_max = new_max
    if cur:
        batches.append(cur)

    scored = {}
    with torch.no_grad():
        for b, batch_idx in enumerate(batches):
            data = [(str(i), sequences[i][:TRUNCATION_SEQ_LENGTH]) for i in batch_idx]
            _, _, toks = batch_converter(data)
            toks = toks.to(device)
            logits = model(toks)["logits"]
            log_probs = _per_position_log_probs(logits)
            for row, i in enumerate(batch_idx):
                L = min(len(sequences[i]), TRUNCATION_SEQ_LENGTH)
                # positions 1..L are residues (0 is BOS, L+1 is EOS)
                token_ids = toks[row, 1 : L + 1]
                lp = log_probs[row, 1 : L + 1, :]
                per_res = lp.gather(1, token_ids.unsqueeze(1)).squeeze(1)
                mean_pll = float(per_res.mean().item())
                scored[i] = {"mean_pll": mean_pll, "n_residues": int(L)}
            print(f"  [swoop] batch {b + 1}/{len(batches)} ({len(batch_idx)} seqs)")

    for i in range(len(sequences)):
        results.append(scored[i])
    return results


def _score_masked(model, alphabet, device, sequences: List[str]) -> List[dict]:
    """Exact masked marginals: mask each position in turn, read log p(true|context).
    O(L) forwards per sequence (positions batched together)."""
    batch_converter = alphabet.get_batch_converter(TRUNCATION_SEQ_LENGTH)
    mask_idx = alphabet.mask_idx
    results: List[dict] = []

    with torch.no_grad():
        for i, seq in enumerate(sequences):
            seq = seq[:TRUNCATION_SEQ_LENGTH]
            L = len(seq)
            _, _, toks = batch_converter([(str(i), seq)])
            toks = toks.to(device)  # (1, L+2)
            true_tokens = toks[0, 1 : L + 1].clone()

            # Build one masked copy per residue position; batch them.
            per_res = torch.empty(L, dtype=torch.float64)
            # Chunk to bound memory on long sequences.
            CHUNK = 64
            pos = 0
            while pos < L:
                chunk_positions = list(range(pos, min(pos + CHUNK, L)))
                batch = toks.repeat(len(chunk_positions), 1)
                for r, p in enumerate(chunk_positions):
                    batch[r, p + 1] = mask_idx
                logits = model(batch)["logits"]
                log_probs = _per_position_log_probs(logits)
                for r, p in enumerate(chunk_positions):
                    per_res[p] = float(log_probs[r, p + 1, true_tokens[p]].item())
                pos += CHUNK
            mean_pll = float(per_res.mean().item())
            results.append({"mean_pll": mean_pll, "n_residues": int(L)})
            print(f"  [masked] {i + 1}/{len(sequences)} ({L} aa)")
    return results


def compute_pseudo_perplexity(
    fasta_path: str,
    *,
    save_path: Optional[str] = None,
    model_location: str = DEFAULT_MODEL,
    method: str = "swoop",
    toks_per_batch: int = 4096,
    nogpu: bool = False,
) -> pd.DataFrame:
    """Score ESM pseudo-perplexity for every sequence in `fasta_path`, writing a
    CSV keyed by ID. method: 'swoop' (fast single-pass approx) or 'masked' (exact)."""
    identifiers, sequences = separate_identifiers(
        load_fasta_sequences(fasta_path, load_identifiers=True)
    )
    identifiers = [i.split(" ", 1)[0] for i in identifiers]
    print(f"Read {fasta_path} with {len(sequences)} sequences")

    model, alphabet = pretrained.load_model_and_alphabet(model_location)
    model.eval()
    device = "cpu"
    if torch.cuda.is_available() and not nogpu:
        model = model.cuda()
        device = "cuda"
        print("Transferred model to GPU")

    if method == "swoop":
        scored = _score_swoop(model, alphabet, device, sequences, toks_per_batch)
    elif method == "masked":
        scored = _score_masked(model, alphabet, device, sequences)
    else:
        raise ValueError(f"Unknown method '{method}' (expected 'swoop' or 'masked').")

    rows = []
    for ident, s in zip(identifiers, scored):
        mean_pll = s["mean_pll"]
        ppl = float(np.exp(-mean_pll)) if np.isfinite(mean_pll) else np.nan
        rows.append(
            {
                "ID": ident,
                "esm_pseudo_perplexity": ppl,
                "esm_mean_pll": mean_pll,
                "n_residues": s["n_residues"],
            }
        )
    df = pd.DataFrame(rows)[COLUMNS]

    if save_path is None:
        partial = os.path.splitext(fasta_path)[0]
        save_path = partial + "_esm_pseudo_perplexity.csv"
    df.to_csv(save_path, index=False)
    print(f"Wrote {len(df)} rows to {save_path}")
    return df
