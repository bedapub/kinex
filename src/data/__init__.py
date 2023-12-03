from importlib import resources
import json
import pandas as pd

def get_pssm() -> pd.DataFrame:
    with resources.path("kinex.data", "pssm_table.csv") as df:
        return pd.read_csv(df, index_col=0)
    
def get_groups() -> dict:
    with resources.path("kinex.data", "groups.json") as file_path:
        with open(file_path) as json_file:
            return json.load(json_file)
        
def get_experiments() -> dict:
    with resources.path("kinex.data", "experiments.json") as file_path:
        with open(file_path) as json_file:
            return json.load(json_file)