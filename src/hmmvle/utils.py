import re

import pandas as pd

from Bio import SeqIO


def process_hyddb(filepath: str) -> pd.DataFrame:

    metadata_df = []

    for seq in SeqIO.parse(seq_file, format="fasta"):

        # Fix missing separator in ID
        seq_id = re.sub(
            pattern="\\.[0-9]+(_)",
            repl="\\.[0-9]+(\\|)",
            string=str(seq.id)
        )

        # Fix trailing "-"
        seq_seq = re.sub(
            pattern="-$",
            repl="",
            string=str(seq.seq)
        )

        seq_id, seq_species, seq_group = seq_id.split("|")
        seq_class = re.findall(r"\[([A-Za-z]+)\]", seq_group)
        seq_group = \
            seq_group.split("_Group_")[-1] \
            if "_Group_" in seq_group else ""

        # Check for more than one occurrence between brackets
        if len(seq_class) != 1: raise NotImplementedError
        else: seq_class = seq_class[0]

        # For Fe class, remove group
        seq_group = f"{seq_class}-{seq_group}" if len(seq_group) else seq_class

        metadata_df.append(
            pd.Series({
                "id": seq_id,
                "species": seq_species,
                "group": seq_group,
                "seq": str(seq_seq)
            }).to_frame().T
        )

    return pd.concat(metadata_df)


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
