import numpy as np
import pandas as pd

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

    msa_tp = build_msa(tp_df)
    msa_fp = build_msa(fp_df)

    # Get IDs for later filtering
    tp_ids = [seq.id.decode() for seq in msa_tp]
    fp_ids = [seq.id.decode() for seq in msa_fp]

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


def hmm_mode(
    hmm_tp: pyhmmer.plan7.HMM,
    hmm_fp: pyhmmer.plan7.HMM,
    alphabet: pyhmmer.easel.Alphabet
) -> pyhmmer.plan7.HMM:
    """
    Algorithm adapted from:

    * Srivastava, P. K., Desai, D. K., Nandi, S., & Lynn, A. M. (2007). 
    HMM-ModE - Improved classification using profile hidden Markov models by 
    optimising the discrimination threshold and modifying emission 
    probabilities with negative training sequences. BMC Bioinformatics, 8(1), 
    104. https://doi.org/10.1186/1471-2105-8-104

    * Sinha, S., & Lynn, A. (2014). HMM-ModE: Implementation, benchmarking and validation with HMMER3. BMC Research Notes, 7(1), 483. https://doi.org/10.1186/1756-0500-7-483
    """

    assert hmm_tp.insert_emissions.shape == hmm_fp.insert_emissions.shape, \
        "[ERROR] Number of nodes must be the same in both TP and FP models!"

    background = pyhmmer.plan7.Background(alphabet)

    n_nodes, n_amino = hmm_tp.match_emissions.shape

    # Avoid first line as per https://pyhmmer.readthedocs.io/en/stable/api/plan7/hmms.html#pyhmmer.plan7.HMM.match_emissions
    for i in range(1, n_nodes):

        # Get relative entropy for both false positives and null
        re_pos_by_neg  = 0.0
        re_pos_by_null = 0.0

        for j in range(n_amino):

            # Avoid division by zero
            fp_denominator = hmm_fp.match_emissions[i, j]
            if fp_denominator == 0.0: fp_denominator = 1e-20

            re_pos_by_neg += hmm_tp.match_emissions[i, j] * \
                np.log2(
                    hmm_tp.match_emissions[i, j] / \
                    fp_denominator
                )

            re_pos_by_null += hmm_tp.match_emissions[i, j] * \
                np.log2(
                    hmm_tp.match_emissions[i, j] / \
                    background.residue_frequencies[j]
                )

        if re_pos_by_neg > re_pos_by_null:

            new_p = np.zeros(n_amino)
            new_null = np.zeros(n_amino)
            denominator = 0

            for k in range(n_amino):

                # Define new null
                # if background.residue_frequencies[k] > hmm_fp.match_emissions[i, k]:
                #     new_null[k] = background.residue_frequencies[k]
                # else:
                #     new_null[k] = hmm_fp.match_emissions[i, k]
                # Same as the following
                new_null[k] = max(
                    background.residue_frequencies[k],
                    hmm_fp.match_emissions[i, k]
                )

                # If P(null) > P(FP)
                if background.residue_frequencies[k] > hmm_fp.match_emissions[i, k]:
                    new_p[k] = hmm_tp.match_emissions[i, k]
                else:
                    new_p[k] = hmm_tp.match_emissions[i, k] * \
                        background.residue_frequencies[k] / \
                        new_null[k]

                denominator += new_p[k]

            # Get final probabilities
            new_score = -np.log2(new_p / denominator)
            final_p = np.exp2(-new_score)

            # Update HMM profile
            for k in range(len(final_p)):
                hmm_tp.match_emissions[i, k] = final_p[k]

            # Validate HMM profile
            hmm_tp.validate(tolerance=1e-4)

    return hmm_tp
