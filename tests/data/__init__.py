from importlib import resources
import pandas as pd

def get_test_pssm() -> pd.DataFrame:
    source = files("tests.data").joinpath('test_pssm_table.csv')
    with as_file(source) as df:
        return pd.read_csv(df, index_col=0)
    
def get_test_scoring_matrix() -> pd.DataFrame:
    source = files("tests.data").joinpath('test_scoring_matrix.csv')
    with as_file(source) as df:
        return pd.read_csv(df, index_col=0)
    
def get_test_input_sites() -> pd.DataFrame:
    source = files("tests.data").joinpath('test_input_sites.csv')
    with as_file(source) as df:
        return pd.read_csv(df)
    

from importlib_resources import files, as_file