import pandas as pd

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
