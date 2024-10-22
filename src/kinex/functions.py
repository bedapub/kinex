import pandas as pd
import numpy as np
from math import sqrt, pow


def get_sequence_format(sequence: str) -> str:
    """
    Returns the sequence format.

    Parameters
    ----------
    sequence : str
        Phosphosite sequence

    Returns
    -------
    str 
        '*'. Phosphorylation site is marked with an asterix "*", on the right side of the phospho-acceptor: GRNS*PLs*PVQA
        '(ph)'. Phosphorylation site is marked with "(ph)", on the right side of the phospho-acceptor: GRNSLs(ph)PVQA
        'central'. Phosphorylation site is in the middle of the sequence; the lenght of the sequence is odd: GRNSLSPVQAS
        'unsupported format'. None of the above-mentioned
    """

    if '*' in sequence:
        return '*'
    elif '(ph)' in sequence:
        return '(ph)'
    else:
        if len(sequence) % 2 == 1:
            return 'central'
    return 'unsupported'


def check_sequence(sequence: str, sequence_format: str) -> bool:
    """
    Checks the validity of the sequence. 

    The input format should be as follows:
    - Phosphorylation site is marked with an asterix "*", on the right side of the phospho-acceptor: GRNSLs*PVQA or,
    - Phosphorylation site is in the middle of the sequence; the lenght of the sequence is odd.
    - Multiple phosphorylation sites are possible. Eg: GRNS*LPs*PVQA.

    - The site can be made up of any of the 20 amino acids, 'X' to hide a position, and '_' to show truncation: GRXNSLs*P_VQA.

    Parameters
    ----------
    sequence : str
        Phosphosite sequence
    sequence_format : str
        Format of the sequence, can be '*' or 'central'

    Returns
    -------
    Boolean 
        True.  If the input string meets the requirements.
        False. If the input string does not meet the requirements.
    """

    allowed_characters = (
        "P",
        "G",
        "A",
        "C",
        "S",
        "T",
        "V",
        "I",
        "L",
        "M",
        "F",
        "Y",
        "W",
        "H",
        "K",
        "R",
        "Q",
        "N",
        "D",
        "E",
        "s",
        "t",
        "y",
        "X",
        "_",
        "*",
    )

    if sequence_format == "*":
        if "**" in sequence:
            return False
        if not any([x in sequence for x in ["S*", "T*", "s*", "t*"]]):
            return False
        for aminoacid in sequence:
            if aminoacid not in allowed_characters:
                return False
    elif sequence_format == "central":
        if sequence[len(sequence) // 2] not in ("S", "T", "s", "t"):
            return False
        for aminoacid in sequence:
            if aminoacid not in [x for x in allowed_characters if x != "*"]:
                return False
    return True


def get_columns(sequence: str, sequence_format: str = "*") -> list:
    """
    Makes a list of column names based on the aminoacid and position in the input sequence. 
    With the phospho-priming option, it includes the phsopho-residues in the phospho-acceptor's vicinity. 

    Parameters
    ----------
    sequence : str
        Phosphosite sequence
    sequence_format : str, default '*'
        Specify which sequence format to use: '*' or 'central'

    Returns
    -------
    list
        List of strings corresponding to the PSSM columns names, based on the position-aminoacid relation in the input sequence.

    Example
    -------
    >>> get_columns(sequence='GRNtSLs*PVQA', sequence_format="*", phospho_priming=True)
    ['-5R', '-4N', '-3t', '-2S', '-1L', '1P', '2V', '3Q', '4A']
    """

    sequence_length = len(sequence)
    column = []
    part_id = 0

    if sequence_format == "central":
        parts = [
            sequence[: sequence_length // 2],
            sequence[sequence_length // 2 + 1:],
        ]  # split the word in half
        for part in parts:  # take first half and second half of the word
            if part_id == 0:
                part = part[::-1]  # remove the S or T from position 0
            for position, aminoacid in enumerate(part):
                if (
                    aminoacid == "_" or aminoacid == "X"
                ):  # jump over a missing character ("_") or trucation ("X").
                    continue
                if part_id == 0:
                    pos = (position + 1) * (-1)
                else:
                    pos = position + 1
                if pos in range(-5, 5):
                    column.append(f"{pos}{aminoacid}")
            part_id += 1
    elif sequence_format == "*":
        parts = sequence.split("*")
        for part in parts:  # take first half and second half of the word
            if part_id == 0:
                part = part[::-1][1:]  # remove the S or T from position 0
            for position, aminoacid in enumerate(part):
                if (
                    aminoacid == "_" or aminoacid == "X"
                ):  # jump over a missing character ("_") or trucation ("X").
                    continue
                if part_id == 0:
                    pos = (position + 1) * (-1)
                else:
                    pos = position + 1
                if pos in range(-5, 5):
                    column.append(f"{pos}{aminoacid}")
            part_id += 1

    return sorted(column, key=lambda item: int(item[:-1]))


def score(sequence: str, sequence_format: str, pssm: pd.DataFrame, favorability: bool = False) -> pd.DataFrame:
    """
    Computes the scores for each of the 303 kinases present in the PSSM table using the list of columns returned by get_columns function.  

    Parameters
    ----------
    sequence : str
        Phosphosite sequence
    sequence_format : str, default '*'
        Specify which sequence format to use: '*' or 'central'
    pssm: pandas.DataFrame
        Position Specific Score Matrix of 303 Ser/Thr kinases
    favorability: bool, default False
        Enable/Disable phospho-acceptor favorability

    Returns
    -------
    pandas.DataFrame
        A table of length 303 with index column 'kinase' and a 'score' column containing the calculated scores.

    Example
    -------
    >>> score(sequence='EGRNSLS*PVQATQ', sequence_format='*', pssm=pssm, phospho_priming=False, favorability=False)
                score
    kinase           
    NLK      8.432319
    JNK3    19.764808
    CDK4    12.719801
    SMG1     1.286346
    JNK1    19.502893
    ...           ...
    TLK2     0.056653
    RAF1     0.132770
    PASK     0.118487
    CDC7     0.160988
    TLK1     0.024576

    [303 rows x 1 columns]
    """
    pssm = pssm.reset_index()

    columns = get_columns(
        sequence, sequence_format)
    columns.append("kinase")

    if favorability == True:
        seq_upper = sequence.upper()
        if seq_upper[len(sequence) // 2] == "S" or "S*" in seq_upper:
            columns.append("0S")
        elif seq_upper[len(sequence) // 2] == "T" or "T*" in seq_upper:
            columns.append("0T")

    df = pssm[columns]
    df.insert(0, "score", df.prod(axis=1, numeric_only=True))
    df = df[["kinase", "score"]]
    df = df.set_index("kinase")
    return df


def get_distances(experiment1, experiment2):
    enrich = np.array(experiment1['dominant_enrichment_value_log2']) - \
        np.array(experiment2['dominant_enrichment_value_log2'])
    p_val = np.array(experiment1['dominant_p_value_log10_abs']) - \
        np.array(experiment2['dominant_p_value_log10_abs'])
    return np.power(np.power(enrich, 2) + np.power(p_val, 2), 0.5)