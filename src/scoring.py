import pandas as pd

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