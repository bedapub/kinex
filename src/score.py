import pandas as pd

class Score:
    """
    Score of the site sequence based on pssm

    Attributes
    ----------
    sequence : str
        A string representing a sequence of Amino-Acids
    scores : pandas.DataFrame
        containing scores, log2(scores) and percentiles for each kinase
    median_percentile : int
        Median value of all percentiles

    Methods
    -------
    promiscuity_index : 
    """

    def __init__(self, sequence: str, ranking: pd.DataFrame) -> None:
        # TODO comment fucntion
        self.sequence = sequence
        # self.columns = columns
        self.ranking = ranking
        self.median_percentile = self.ranking["percentile_score"].median()
        # Calculate the median percentile and promiscuity_index

    def __repr__ (self):
        return f"Scoring results for {self.sequence}"
    
    
    # TODO check attributes format before setting
    @property
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, sequence):
        self._sequence = sequence

    @property
    def columns(self):
        return self._columns
    @columns.setter
    def columns(self, columns):
        self._columns = columns

    @property
    def ranking(self):
        return self._ranking
    @ranking.setter
    def ranking(self, ranking):
        self._ranking = ranking

    def promiscuity_index(self, limit: int = 90) -> int:
        return self.ranking["percentile_score"][self.ranking["percentile_score"] > limit].count()
    
    def top(self, number: int = 15):
        return self.ranking.head(number)