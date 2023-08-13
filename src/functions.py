import pandas as pd

def get_sequence_format(sequence: str) -> str:
    """
    Returns the sequence format

    Parameters
    ----------
    sequence : str
        Phosphosite sequence

    Returns
    -------
    str 
        '*'. Phosphorylation site is marked with an asterix "*", on the right side of the phospho-acceptor: GRNSLs*PVQA
        'central'. Phosphorylation site is in the middle of the sequence, the lenght of the sequence is odd 
        'unsupported format'. None of above mentioned
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
    Checks if the input file contains the right sequence format

    The input format should be as follows:
    - Phosphorylation site is marked with an asterix "*", on the right side of the phospho-acceptor: GRNSLs*PVQA
    - The site sequence can include any of the 20 amino acids, 'X' for masking a position, and '_' for truncation: GRXNSLs*P_VQA
    - To include phosphorylated residues (pS/pT/pY), use lower-case letters (s/t/y) and 
    PSVEPPLs*QEtFSDL (with phospho-priming option checked)

    Parameters
    ----------
    sequence : str
        Phosphosite sequence
    sequence_format : str
        Format of the sequence, can be '*' or 'central'

    Returns
    -------
    Boolean 
        True.  If the input string meets the requirements
        False. If the input string does not meet the requirements 

    Raises
    ------
    ValueError
        If sequence format is not supported
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

def get_columns(sequence: str, sequence_format: str = "*", phospho_priming: bool = False) -> list:
    """
    Makes a list of column names required for the calculation depending on the inputted sequence.

    Parameters
    ----------
    sequence : str
        Phosphosite sequence
    sequence_format : str, default '*'
        To tell which format to use: '*' or 'central'

    Returns
    -------
    list
        List of column names of the pssm table required for the calculation, depending on the inputted sequence

    Raises
    ------
    ValueError
        If sequence_format not supported
    """

    sequence_length = len(sequence)
    column = []
    part_id = 0

    if not phospho_priming:
        sequence = sequence.upper()

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

def score(sequence: str, sequence_format: str, pssm: pd.DataFrame, phospho_priming: bool = False) -> pd.DataFrame:
    columns = get_columns(
            sequence, sequence_format, phospho_priming)
    columns.append("kinase")
    df = pssm[columns]
    df.insert(0, "score", df.prod(axis=1, numeric_only=True))
    df = df[["kinase", "score"]]
    df = df.set_index("kinase")
    return df