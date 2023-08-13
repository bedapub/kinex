from importlib import resources
import pandas as pd

def get_pssm() -> pd.DataFrame:
    with resources.path("data", "pssm_table.csv") as df:
        return pd.read_csv(df, index_col=0)