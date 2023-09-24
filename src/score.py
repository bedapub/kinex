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

    def __init__(self,sequence: list, ranking: list) -> None:
        self.sequence = sequence
        self.n_sequences = len(ranking)
        median_percentile = []
        if self.n_sequences == 1:
            self.ranking = ranking[0]
            self.median_percentile = ranking[0]["percentile_score"].median()
        else:
            self.ranking = ranking
            for score in ranking:
                median_percentile.append(score["percentile_score"].median())
            self.median_percentile = median_percentile

    def __repr__ (self):
        return f"Scoring results for {self.sequence}"

    def promiscuity_index(self, limit: int = 90, subsequence: int = 0) -> int:
        if self.n_sequences == 1:
            return self.ranking["percentile_score"][self.ranking["percentile_score"] > limit].count()
        else:
            return self.ranking[subsequence]["percentile_score"][self.ranking[subsequence]["percentile_score"] > limit].count()
    
    def top(self, number: int = 15, subsequence: int = 0):
        if self.n_sequences == 1:
            return self.ranking.head(number)
        else:
            return self.ranking[subsequence].head(number)