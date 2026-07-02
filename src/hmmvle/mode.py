import pyhmmer
import pyfamsa
import pytrimal

from hmmvle.hmm import build_msa


def build_profiles_aligned(
    tp_df: pd.DataFrame,
    fp_df: pd.DataFrame,
    name: str,
    alphabet: pyhmmer.easel.Alphabet
) -> list[pyhmmer.plan7.HMM]:

    background = pyhmmer.plan7.Background(alphabet)

    msa_tp = _build_msa(tp_df)
    msa_fp = _build_msa(fp_df)

    # Get IDs for later filtering
    tp_ids = [seq.id for seq in msa_tp]
    fp_ids = [seq.id for seq in msa_fp]

    # Align both TP and FP MSAs
    aligner = pyfamsa.Aligner(
        guide_tree="nj",
        keep_duplicates=False
    )
    msa = aligner.align_profiles(
        profile1=msa_tp,
        profile2=msa_fp
    )

    # ------------------------------------------------------------------------ #
    # Trimming is required to guarantee same number of nodes in both HMMs

    alignment = pytrimal.Alignment(
        names=[seq.id for seq in msa],
        sequences=[seq.sequence for seq in msa]
    )
    trimmer = pytrimal.AutomaticTrimmer(method="gappyout")

    trimmed = trimmer.trim(alignment)

    trimmed_msa = trimmed.to_pyhmmer()
    trimmed_msa.name = "combined".encode()

    # ------------------------------------------------------------------------ #

    msa_text_tp  = pyhmmer.easel.TextMSA(
        name=f"{name}".encode(),
        sequences=[
            pyhmmer.easel.TextSequence(
                name=name,
                sequence=seq
            )
            for name, seq in zip(trimmed_msa.names, trimmed_msa.alignment)
            if name in tp_ids
        ]
    )
    msa_text_fp  = pyhmmer.easel.TextMSA(
        name=f"{name}-FP".encode(),
        sequences=[
            pyhmmer.easel.TextSequence(
                name=name,
                sequence=seq
            )
            for name, seq in zip(trimmed_msa.names, trimmed_msa.alignment)
            if name in fp_ids
        ]
    )

    # WARNING: set symfrac to 0.0 so that all columns in both MSA are assigned 
    # as consensus. Otherwise, mismatches can happen when different consensus 
    # columns are assigned to different MSAs.
    builder = pyhmmer.plan7.Builder(
        architecture="fast",
        alphabet=alphabet,
        symfrac=0.0,
        # fragthresh=0.0, # Set to 0.0. if none of the sequences are fragments
    )

    hmm_tp, _, _ = builder.build_msa(msa_text_tp.digitize(alphabet), background)
    hmm_fp, _, _ = builder.build_msa(msa_text_fp.digitize(alphabet), background)

    return hmm_tp, hmm_fp
