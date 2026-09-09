import re
from io import StringIO

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


def get_hits(
    filepath: str,
    sep: str = "$"
) -> pd.DataFrame:

    hmmer_colnames = [
        "target_name",
        "target_accession",
        "query_name",
        "query_accession",
        "e_value_full_seq",
        "score_full_seq",
        "bias_full_seq",
        "e_value_best_dom",
        "score_best_dom",
        "bias_best_dom",
        "exp",
        "reg",
        "clu",
        "ov",
        "env",
        "dom",
        "rep",
        "inc",
        "description_of_target"
    ]

    with open(filepath, mode="r") as handle:
        hits = handle.readlines()

        # Remove comment lines
        hits = [line for line in hits if not line.startswith("#")]

        # Remove whitespaces in description of target column
        hits = [re.sub("#\\s.*\\sID", "ID", line) for line in hits]

        # Replace multiple whitespaces by separator
        hits = [re.sub("\\s+", sep, line) for line in hits]

        # Manage potential whitespaces in description
        hits = [
            "$".join(
                line.split(sep)[:len(hmmer_colnames) - 1] + \
                ["_".join(line.split(sep)[len(hmmer_colnames) - 1:])]
            )
            for line in hits
        ]

        # Remove trailing separator
        hits = [line.rstrip(sep) for line in hits]

        # Concatenate by new line character
        hits = "\n".join(hits)

    return pd.read_csv(
        StringIO(hits),
        names=hmmer_colnames,
        sep=sep
    )
