import pandas as pd

class Score:
    """
    Score of the site sequence based on pssm

    Parameters
    ----------
    sequence : str
        A string representing a sequence of Amino-Acids
    columns : list
        List of column names for the sequence of the pssm table required for the calculation
    scores : pandas.DataFrame
        containing scores, log2(scores) and percentiles for each kinase
    median_percentile : int
        Median value of all percentiles

    Methods
    -------
    promiscuity_index : 
    """

    def __init__(self, sequence: str, scores: pd.DataFrame) -> None:
        # TODO comment fucntion
        self.sequence = sequence
        # self.columns = columns
        self.scores = scores
        self.median_percentile = self.scores["percentile_score"].median()
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
    def scores(self):
        return self._scores
    @scores.setter
    def scores(self, scores):
        self._scores = scores

    def promiscuity_index(self, limit: int = 90) -> int:
        return self.scores["percentile_score"][self.scores["percentile_score"] > limit].count()
    
    def top(self, number: int = 15):
        return self.scores.head(number)