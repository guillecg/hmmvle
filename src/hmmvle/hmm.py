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


def get_seqs(filepath: str) -> pd.DataFrame:

    seqs = []

    with open(filepath, mode="r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs.append(
                pd.Series({
                    "seq_id": record.id,
                    "seq": "".join(record.seq)
                }).to_frame().T
            )

    return pd.concat(seqs)


def define_hits(
    row: pd.Series,
    hmm_thr: float,
    group: str
) -> str:
    # There are no false negatives because by definition the threshold includes
    # all the considered true positives (threshold is minimum among them)
    if row["score_full_seq"] >= hmm_thr and row["group"] == group:
        return "True positive"
    elif row["score_full_seq"] >= hmm_thr and row["group"] != group:
        return "False positive"
    else:
        return "True negative"
