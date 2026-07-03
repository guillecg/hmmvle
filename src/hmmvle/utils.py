import pandas as pd


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
