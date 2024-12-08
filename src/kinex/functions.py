import requests
from importlib import resources
from pathlib import Path

import numpy as np

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

def get_distances(experiment1: dict, experiment2: dict) -> np.ndarray:
    """
    
    Calculate the Euclidean distance between two experiments based on their
    dominant enrichment and p-value scores.
    
    """
    
    # Validate input keys
    required_keys = ['dominant_enrichment_value_log2', 'dominant_p_value_log10_abs']
    for key in required_keys:
        if key not in experiment1 or key not in experiment2:
            raise ValueError(f"Both experiments must contain the key '{key}'.")

    # Convert input values to NumPy arrays
    enrich1 = np.array(experiment1['dominant_enrichment_value_log2'])
    enrich2 = np.array(experiment2['dominant_enrichment_value_log2'])
    p_val1 = np.array(experiment1['dominant_p_value_log10_abs'])
    p_val2 = np.array(experiment2['dominant_p_value_log10_abs'])

    # Ensure arrays have the same length
    if enrich1.shape != enrich2.shape or p_val1.shape != p_val2.shape:
        raise ValueError("The arrays in the two experiments must have the same shape.")

    # Calculate the Euclidean distance
    enrich_diff = enrich1 - enrich2
    p_val_diff = p_val1 - p_val2
    distances = np.sqrt(np.power(enrich_diff, 2) + np.power(p_val_diff, 2))
    
    return distances


def download_file_to_resource(url: str, resource_name: str) -> None:
    """
    Downloads a file from a given URL and saves it to the specified resource path.

    Parameters:
    ----------
    url : str
        The URL of the file to be downloaded.
    resource_name : str
        The name of the resource file to save (e.g., "default_scoring_matrix_tyr.csv.gz").

    Raises:
    -------
    ValueError:
        If the URL or resource details are invalid.
    requests.exceptions.RequestException:
        If there are issues with the HTTP request (e.g., network error, 404).
    IOError:
        If there's an error saving the file.
    """
    try:     
        # Determine the file save path using importlib.resources
        with resources.path("kinex.resources", resource_name) as file_path:
            save_path = Path(file_path)

        print(f"Starting download from: {url}")
        response = requests.get(url, stream=True, timeout=10)
        
        # Raise an error for HTTP codes 4xx/5xx
        response.raise_for_status()
        
        
        with open(save_path, "wb") as file:
            file.write(response.content)
        
        print(f"File successfully downloaded and saved to: {save_path}")
    
    except requests.exceptions.MissingSchema:
        raise ValueError("The provided URL is not valid.")
    except requests.exceptions.RequestException as e:
        print(f"Error during file download: {e}")
        raise
    except IOError as e:
        print(f"Error saving the file to {save_path}: {e}")
        raise
