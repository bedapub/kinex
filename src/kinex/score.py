from kinex.sequence import Sequence


class Score:
    """
    Score of the site validation based on pssm

    Attributes
    ----------
    sequence : str
        A string representing a validation of aminoacids
    ranking : pandas.DataFrame
        containing scores, log2(scores) and percentiles for each kinase
    median_percentile : int
        Median value of all percentile scores for a validation
        
    Methods
    -------
    promiscuity_index(self, limit: int = 90) -> int
        The number of kinases scoring above the {limit}th percentile

    top(self, number: int = 15) -> 
        Top {number} ranked kinases according to percentile score. 

    """

    def __init__(self, sequence: Sequence, ranking: list) -> None:
        self.sequence = sequence
        self.ranking = ranking

        median_percentile = []
        for score in ranking:
            median_percentile.append(score["percentile_score"].median())
        self.median_percentile = median_percentile

    def __repr__(self):
        return f"Scoring results for {self.sequence.sequence_string}"

    def promiscuity_index(self, limit: int = 90) -> list:
        promiscuity_index = []
        for score in self.ranking:
            promiscuity_index.append(score["percentile_score"][score["percentile_score"] > limit].count())
        return promiscuity_index

    def top(self, number: int = 15) -> list:
        top = []
        for score in self.ranking:
            top.append(score.head(number))
        return top