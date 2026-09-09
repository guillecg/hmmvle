import pandas as pd

import pyhmmer


def search_single(
    seq_path: str,
    hmm: pyhmmer.plan7.HMM,
    alphabet: pyhmmer.easel.Alphabet
) -> pd.DataFrame:

    background = pyhmmer.plan7.Background(alphabet)
    pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)

    with pyhmmer.easel.SequenceFile(
        seq_path,
        digital=True,
        alphabet=alphabet
    ) as handle:
        hits = pipeline.search_hmm(
            query=hmm,
            sequences=handle
        )

    return hits
