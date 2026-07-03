import re


def _fix_missing_separator(seq_id: str) -> str:
    return re.sub(
        pattern=r"(\.[0-9]+)_",
        repl=r"\1|",
        string=str(seq_id)
    )


def _fix_trailing_dash(seq: str) -> str:
    return re.sub(
        pattern="-$",
        repl="",
        string=str(seq)
    )
