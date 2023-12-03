import pandas as pd
from importlib_resources import files, as_file

def get_test_pssm() -> pd.DataFrame:
    source = files("kinex.tests.data").joinpath('test_pssm_table.csv')
    with as_file(source) as df:
        return pd.read_csv(df, index_col=0)
    
def get_test_scoring_matrix() -> pd.DataFrame:
    source = files("kinex.tests.data").joinpath('test_scoring_matrix.csv')
    with as_file(source) as df:
        return pd.read_csv(df, index_col=0)
    
def get_test_input_sites() -> pd.DataFrame:
    source = files("kinex.tests.data").joinpath('test_input_sites.csv')
    with as_file(source) as df:
        return pd.read_csv(df)

