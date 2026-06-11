"""Extract EnzymeExplorer (EE) ESM-1v-TPS PLM embeddings — the PLM feature block of
EE's production model ``PlmDomainsRandomForest__tps_esm-1v-subseq_..._domains_subset``.

This is **block (B)** of the two EE feature blocks (the structure/function domain block
is produced by the sibling ``extract_ee_domain_features.py``). It reproduces the exact
PLM representation that EE's ``scripts/easy_predict.py`` feeds to the production
classifier:

  * Model: ESM-1v (``esm1v_t33_650M_UR90S_1``) with the TPS-finetuned "subseq"
    checkpoint ``checkpoint-tps-esm1v-t33-subseq.ckpt`` (EE's
    ``esm_transformer_utils.get_model_and_tokenizer("esm-1v-finetuned-subseq")``).
  * Representation layer 33, **mean-pooled** over residue tokens (excluding BOS/EOS) ->
    a 1280-d vector per sequence (EE's ``compute_embeddings``, ``max_len=1022``).

The production classifier's ``plm_feat_indices_subset`` is the full ``range(1280)``
(identity) for every fold, so the raw 1280-d mean-pooled embedding IS the PLM block fed
to the model — no column subsetting needed.

Output: one row per protein, columns = ``id`` + ``emb_0..emb_1279``, keyed by
``Enzyme_marts_ID``.

Must run in the EE ``enzyme_explorer`` env (esm + torch). GPU strongly recommended.
Run from the EE ``scripts/`` dir so ``data/plm_checkpoints/...`` resolves to the
production ``scripts/data/`` bundle. Sequences are read from the input CSV directly
(NOT from the structures), matching easy_predict.py for proteins routed through the
domain branch — easy_predict reads the seq from the PDB, but the PDB seq == the input
sequence for these ESMFold structures; we read the canonical input sequence for fidelity
and to cover proteins without a structure.
"""

from __future__ import annotations

import argparse
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

from enzymeexplorer.src.embeddings_extraction.esm_transformer_utils import (
    compute_embeddings,
    get_model_and_tokenizer,
)

EMB_DIM = 1280
PLM_MAX_SEQ_LEN = 1022  # easy_predict.py default --plm-max-seq-len


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract EE ESM-1v-TPS-subseq mean-pooled PLM embeddings (1280-d) "
        "for a set of protein sequences."
    )
    parser.add_argument("--sequences_csv", required=True,
                        help="CSV with --id_column and --sequence_column.")
    parser.add_argument("--id_column", default="Enzyme_marts_ID")
    parser.add_argument("--sequence_column", default="Aminoacid_sequence")
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--max_seq_len", type=int, default=PLM_MAX_SEQ_LEN)
    args = parser.parse_args()

    df_in = pd.read_csv(args.sequences_csv)
    ids = df_in[args.id_column].astype(str).tolist()
    seqs = df_in[args.sequence_column].astype(str).tolist()
    print(f"{len(ids)} sequence(s) from {args.sequences_csv}", flush=True)

    model, batch_converter, alphabet = get_model_and_tokenizer(
        "esm-1v-finetuned-subseq", return_alphabet=True
    )
    embed = partial(
        compute_embeddings,
        bert_model=model,
        converter=batch_converter,
        padding_idx=alphabet.padding_idx,
        model_repr_layer=33,
        max_len=args.max_seq_len,
    )

    all_emb = np.zeros((len(ids), EMB_DIM), dtype=np.float32)
    for start in range(0, len(ids), args.batch_size):
        batch_seqs = seqs[start : start + args.batch_size]
        # truncate exactly as easy_predict.py does for the domain branch
        batch_seqs = [s[: (args.max_seq_len - 2)] if len(s) > args.max_seq_len else s
                      for s in batch_seqs]
        enc, _ = embed(input_seqs=batch_seqs)
        all_emb[start : start + len(batch_seqs)] = enc
        print(f"  embedded {min(start + len(batch_seqs), len(ids))}/{len(ids)}",
              flush=True)

    cols = [f"emb_{i}" for i in range(EMB_DIM)]
    out = pd.DataFrame(all_emb, columns=cols)
    out.insert(0, "id", ids)
    Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(out)} rows x {EMB_DIM} PLM dims to {args.output_csv}", flush=True)


if __name__ == "__main__":
    main()
