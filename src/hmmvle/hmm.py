import re
from io import StringIO

import pandas as pd

from Bio import SeqIO

import pyhmmer
import pyfamsa


def build_msa(df: pd.DataFrame) -> pyfamsa.Alignment:
    sequences = [
        pyfamsa.Sequence(row["id"].encode(), row["seq"].encode())
        for _, row in df.iterrows()
    ]
    aligner = pyfamsa.Aligner(
        guide_tree="nj",
        keep_duplicates=False
    )
    msa = aligner.align(sequences)

    return msa


def build_profile(
    df: pd.DataFrame,
    name: str,
    alphabet: pyhmmer.easel.Alphabet
) -> pyhmmer.plan7.HMM:

    msa = build_msa(df)

    # Covert to TextMSA to digitize it
    msa_text  = pyhmmer.easel.TextMSA(
        name=name.encode(),
        sequences=[
            pyhmmer.easel.TextSequence(
                name=seq.id,
                sequence=seq.sequence.decode()
            )
            for seq in msa
        ]
    )

    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)

    hmm, _, _ = builder.build_msa(msa_text.digitize(alphabet), background)

    return hmm
