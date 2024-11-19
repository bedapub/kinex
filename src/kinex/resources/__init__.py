from importlib import resources
import json
import pandas as pd

def get_pssm_ser_thr() -> pd.DataFrame:
    with resources.path("kinex.resources", "pssm_table_ser_thr.csv") as df:
        return pd.read_csv(df, index_col=0)


def get_pssm_tyr() -> pd.DataFrame:
    with resources.path("kinex.resources", "pssm_table_tyr.csv") as df:
        return pd.read_csv(df, index_col=0)


def get_ser_thr_family() -> dict:
    with resources.path("kinex.resources", "ser_thr_family.json") as file_path:
        with open(file_path) as json_file:
            return json.load(json_file)


def get_ser_thr_family_colors() -> dict:
    with resources.path("kinex.resources", "ser_thr_family_colors.json") as file_path:
        with open(file_path) as json_file:
            return json.load(json_file)


def get_tyr_family() -> dict:
    with resources.path("kinex.resources", "tyr_family.json") as file_path:
        with open(file_path) as json_file:
            return json.load(json_file)


def get_tyr_family_colors() -> dict:
    with resources.path("kinex.resources", "tyr_family_colors.json") as file_path:
        with open(file_path) as json_file:
            return json.load(json_file)


def get_experiments() -> dict:
    with resources.path("kinex.resources", "experiments.json") as file_path:
        with open(file_path) as json_file:
            return json.load(json_file)
        
        
def get_scoring_matrix_ser_thr() -> pd.DataFrame:
    try:
        with resources.path("kinex.resources", "default_scoring_matrix_ser_thr.csv.gz") as file_path:
            return pd.read_csv(file_path, compression='gzip')
    except FileNotFoundError:
         return None

def get_scoring_matrix_tyr() -> pd.DataFrame:
    try:
        with resources.path("kinex.resources", "default_scoring_matrix_tyr.csv.gz") as file_path:
            return pd.read_csv(file_path, compression='gzip')
    except FileNotFoundError:
        return None
