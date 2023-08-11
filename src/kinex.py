import pandas as pd
import numpy as np
import scipy as sc

from input import check_sequence, get_sequence_format
from scoring import get_score

from score import Score
from enrichment import Enrichment

class Kinex:
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

    def score(self, sequence: str, phospho_priming: bool = False, favorability: bool = False, method: str = 'avg'):
        # TODO comment fucntion
        sequence_format = get_sequence_format(sequence)

        # Check method and format
        if not method in ['min', 'max', 'avg']:
            raise ValueError(
                f"Method {method} is not supported. Supported methods: 'min', 'max', 'avg'")
        
        if sequence_format == 'unsupported':
            raise ValueError(f"Sequence format is not supported. Supported formats: '*' and 'central'")

        # Empty list for possible sequences
        sequences = []
        if sequence_format == '(ph)' or sequence_format == '*':
            # Split the sequence to sub-sequences and count them
            sequence = sequence.split(sequence_format)
            num = len(sequence)
            sequence_format = '*'
            
            # Make the last aminoacid lowercase
            for id in range(num):
                if not id == num - 1:
                    sequence[id] = sequence[id][:-1].upper() + sequence[id][-1:].lower()

            # Make possible sequences 
            for id in range(num - 1):
                seq = ''
                for i in range(num):
                    if i == id + 1:
                        seq += "*"
                    seq += sequence[i]
                # Check the validity of subsequence. If invalid, go to next
                if not check_sequence(seq, sequence_format):
                    continue
                sequences.append(seq)

            # If all the subsequences invalid -> sequence invalid
            if len(sequences) == 0:
                raise ValueError("Invalid sequence")
            sequence = sequences

        # Empty dataframe to store scores
        df = pd.DataFrame()
        # Iterate through every sequence for asterisk format
        if sequence_format == "*":
            number_of_seq = len(sequence)
            for id in range(number_of_seq):
                if id == 0:
                    df = get_score(sequence[id], sequence_format, self.pssm, phospho_priming) / number_of_seq if method == 'avg' else get_score(sequence[id], sequence_format, self.pssm, phospho_priming)
                else:
                    if method == 'avg':
                        df = df.add(get_score(sequence[id], sequence_format, self.pssm, phospho_priming) / number_of_seq)
                    else:
                        df = pd.concat([df,get_score(sequence[id], sequence_format, self.pssm, phospho_priming)])
                        if method == 'min':
                            df = df.groupby(df.index).min()
                        elif method == 'max':
                            df = df.groupby(df.index).max()
        # Central format
        elif sequence_format == 'central':
            if not check_sequence(sequence, sequence_format):
                raise ValueError("Invalid sequence")
            df = get_score(sequence, sequence_format, self.pssm, phospho_priming)

        # TODO Add favorability support
        df.insert(1, "log_score", np.log2(df["score"]))

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
        return Score(sequence, df)
    

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
                top15_kinases = self.score(str(input_sites.iloc[id, 0])).top(15).index
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