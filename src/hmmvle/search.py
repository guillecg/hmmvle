import pandas as pd

import pyhmmer
from pyhmmer.easel import SequenceFile


def search_single(
    seq_path: str,
    hmm: pyhmmer.plan7.HMM,
    alphabet: pyhmmer.easel.Alphabet
) -> pd.DataFrame:

    background = pyhmmer.plan7.Background(alphabet)
    pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)
    
    with SequenceFile(seq_path, digital=True, alphabet=alphabet) as seq_file:
        hits = pipeline.search_hmm(hmm, seq_file)

    return hits
