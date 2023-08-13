from importlib import resources
import pandas as pd

def get_test_pssm() -> pd.DataFrame:
    with resources.path("tests.data", "test_pssm_table.csv") as df:
        return pd.read_csv(df, index_col=0)
    
def get_test_scoring_matrix() -> pd.DataFrame:
    with resources.path("tests.data", "test_scoring_matrix.csv") as df:
        return pd.read_csv(df, index_col=0)
    
def get_test_input_sites() -> pd.DataFrame:
    with resources.path("tests.data", "test_input_sites.csv") as df:
        return pd.read_csv(df)