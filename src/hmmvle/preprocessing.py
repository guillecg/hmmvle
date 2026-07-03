import re

from Bio import SeqIO


def _fix_missing_separator(seq_id: str) -> str:
    return re.sub(
        pattern=r"(\.[0-9]+)_",
        repl=r"\1|",
        string=str(seq_id)
    )


def process_hyddb(filepath: str) -> pd.DataFrame:

    metadata_df = []

    for seq in SeqIO.parse(filepath, format="fasta"):

        # Fix missing separator in ID
        seq_id = _fix_missing_separator(seq.id)

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
