import pandas as pd

class Score:
    """
    Score of the site sequence based on pssm

    Attributes
    ----------
    sequence : str
        A string representing a sequence of aminoacids
    ranking : pandas.DataFrame
        containing scores, log2(scores) and percentiles for each kinase
    median_percentile : int
        Median value of all percentile scores for a sequence
        
    Methods
    -------
    promiscuity_index(self, limit: int = 90) -> int
        The number of kinases scoring above the {limit}th percentile

    top(self, number: int = 15) -> 
        Top {number} ranked kinases according to percentile score. 

    """

    def __init__(self, sequence: str, ranking: pd.DataFrame) -> None:
        self.sequence = sequence
        self.ranking = ranking
        self.median_percentile = self.ranking["percentile_score"].median()

    def __repr__ (self):
        return f"Scoring results for {self.sequence}"
    
    @property
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, sequence):
        self._sequence = sequence

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