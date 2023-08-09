import pandas as pd
import numpy as np
import scipy as sc
from statsmodels.stats.multitest import multipletests

from input_check import check_sequence, get_columns, get_sequence_format
from score import Score
from enrichment import Enrichment

class PSSM:
    """
    The class representing a pssm table including the table and methods like scoring and enrichment

    Parameters
    ----------
    pssm : pandas.DataFrame
        Info
    scoring_matrix : pandas.DataFrame
        Info

    Methods
    -------
    get_score : 
    enrichment : 
    """

    def __init__(self, pssm: pd.DataFrame, scoring_matrix: pd.DataFrame):
        """
        Initializes the instance of the PSSM table
        Parameters
        ----------
        pssm : pandas.DataFrame
            A PSSM table containing scores for each position and Amino Acid
        """
        self.pssm = pssm
        self.scoring_matrix = scoring_matrix
        # TODO Check the table format
        # TODO Add favorability attribute

    @property
    def pssm(self):
        return self._pssm

    @pssm.setter
    def pssm(self, pssm):
        self._pssm = pssm

    @property
    def scoring_matrix(self):
        return self._scoring_matrix

    @scoring_matrix.setter
    def scoring_matrix(self, scoring_matrix):
        self._scoring_matrix = scoring_matrix

    def get_score(self, sequence: str, phospho_priming: bool = False, favorability: bool = False):
        # TODO comment fucntion
        # Sequence check
        sequence_format = get_sequence_format(sequence)
        if not check_sequence(sequence):
            raise ValueError("Invalid sequence")

        # Prepare columns
        columns = get_columns(
            sequence, sequence_format,  phospho_priming)
        columns.append("kinase")
        df = self.pssm[columns]
        df.insert(0, "score", df.prod(axis=1, numeric_only=True))

        # TODO Add favorability support
        df.insert(1, "log_score", np.log2(df["score"]))

        df = df[["kinase", "score", "log_score"]]
        df = df.set_index("kinase")

        # Compute percentiles
        percentiles = []
        for kinase in df.index:
            # TODO QUESTION Should we round to 3 decimals?
            percentile = round(
                sc.stats.percentileofscore(
                    self.scoring_matrix[kinase], df.log_score[kinase]
                ),
                3,
            )
            percentiles.append(percentile)
        df.insert(2, "percentile_score", percentiles)
        df = df.sort_values("percentile_score", ascending=False)
        # TODO create an instance of the Score class with attributes sequence, columns, scores
        # return that instance
        return Score(sequence, columns, df)
    

    def enrichment(self, input_sites: pd.DataFrame, fc_threshold: float = 1.5):
        # TODO check input_sites format, make sure there are all necessary columns

        # Empty DataFrame to store the output
        enrichment_table = pd.DataFrame(
            columns=['kinase', 'upregulated', 'downregulated', 'unregulated'])

        # Count the number of sets
        total_upregulated = total_downregulated = total_unregulated = 0
        regulation_list = []
        failed_sites = []


        
        for id in range(len(input_sites)):
            # Get top 15 kinases, check if site is valid
            try:
                top15_kinases = self.get_score(str(input_sites.iloc[id, 0])).top(15).index
            except ValueError:
                failed_sites.append(input_sites.iloc[id, 0])
                regulation_list.append('failed')
                continue
            
            regulation = ""
            
            # TODO check data type conversions
            if float(str(input_sites.iloc[id, 1])) > fc_threshold:
                regulation = "upregulated"
                total_upregulated += 1
            elif float(str(input_sites.iloc[id, 1])) < -fc_threshold:
                regulation = "downregulated"
                total_downregulated += 1
            elif float(str(input_sites.iloc[id, 1])) <= fc_threshold:
                regulation = "unregulated"
                total_unregulated += 1

            regulation_list.append(regulation)

            enrichment_table = pd.concat([enrichment_table, pd.DataFrame(
                {"kinase": top15_kinases, regulation: np.ones(15)})]).groupby('kinase').sum().reset_index()

        # upregulated_adjusted_p_value = multipletests(self.enrichment_table["upregulated_p_value"], method="fdr_bh")
        # self.enrichment_table.insert(6, "upregulated_adjusted_p_value", upregulated_adjusted_p_value[1])

        # Add regulation column to input_sites table
        input_sites.insert(2, 'regulation', regulation_list)
        # TODO implement the Enrichment class and return instance of that class

        # Background adjustment
        if total_unregulated == 0:
            total_unregulated = np.min(
                [total_upregulated, total_downregulated])/2

        return Enrichment(enrichment_table, input_sites, failed_sites, total_upregulated, total_downregulated, total_unregulated)