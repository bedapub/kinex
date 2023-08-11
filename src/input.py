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