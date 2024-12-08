from importlib import resources
import json
import pandas as pd

def get_pssm_ser_thr() -> pd.DataFrame:
    with resources.files("kinex.resources").joinpath("pssm_table_ser_thr.csv").open() as df:
        return pd.read_csv(df, index_col=0)


def get_pssm_tyr() -> pd.DataFrame:
    with resources.files("kinex.resources").joinpath("pssm_table_tyr.csv").open() as df:
        return pd.read_csv(df, index_col=0)


def get_ser_thr_family() -> dict:
    with resources.files("kinex.resources").joinpath("ser_thr_family.json").open() as json_file:
        return json.load(json_file)


def get_ser_thr_family_colors() -> dict:
    with resources.files("kinex.resources").joinpath("ser_thr_family_colors.json").open() as json_file:
        return json.load(json_file)


def get_tyr_family() -> dict:
    with resources.files("kinex.resources").joinpath("tyr_family.json").open() as json_file:
        return json.load(json_file)


def get_tyr_family_colors() -> dict:
    with resources.files("kinex.resources").joinpath("tyr_family_colors.json").open() as json_file:
        return json.load(json_file)


def get_experiments() -> dict:
    with resources.files("kinex.resources").joinpath("experiments.json").open() as json_file:
        return json.load(json_file)
        
        
def get_scoring_matrix_ser_thr() -> pd.DataFrame:
    try:
        with resources.files("kinex.resources").joinpath("default_scoring_matrix_ser_thr.csv.gz").open('rb') as file_path:
            return pd.read_csv(file_path, compression='gzip')
    except FileNotFoundError:
         return None

def get_scoring_matrix_tyr() -> pd.DataFrame:
    try:
        with resources.files("kinex.resources").joinpath("default_scoring_matrix_tyr.csv.gz").open('rb') as file_path:
            return pd.read_csv(file_path, compression='gzip')
    except FileNotFoundError:
        return None

def get_configuration_file() -> dict:
    with resources.files("kinex.resources").joinpath("config.json").open() as json_file:
        return json.load(json_file)
