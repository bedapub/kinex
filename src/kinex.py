import pandas as pd
import numpy as np
# import time
import bisect

from kinex.functions import check_sequence, get_sequence_format, score

from kinex.data import get_pssm

from kinex.score import Score
from kinex.enrichment import Enrichment
from kinex.comparison import Comparison

# import logging
# logging.basicConfig(
#     filename="kinex.log",
#     level=logging.ERROR,
#     format="%(asctime)s %(levelname)s %(message)s",
#     datefmt="%Y-%m-%d %H:%M:%S",
# )


class Kinex:
    """
    The class representing a PSSM table and a scoring matrix needed for scoring and enrichment analysis.

    Attributes
    ----------
    pssm : pandas.DataFrame
        Normalised and scaled densiometries from PSPA experiments. 
        The table cotains on rows the kinases and on columns the positions for each aminoacid.
    scoring_matrix : pandas.DataFrame
        Table containing 82,755 experimentally identified Ser/Thr phosphosites that have been scored by 303 Ser or Thr kinase PSSM.
        The table allows the ranking of kinases, as well as the calculation of promiscuity index and median percentile for each input sequence.

    Methods
    -------
    get_score(self, sequence: str, phospho_priming: bool = False, favorability: bool = False, method: str = 'avg') -> score.Score
        Checks the sequence format and it's validity; computes the scores and ranks the kinases.
    get_enrichment(self, input_sites: pd.DataFrame, fc_threshold: float = 1.5) -> enrichment.Enrichment
        Checks regulation and plots the the enrichment vs p-value. 

    Examples
    --------
    PSSM table
    >>> kinex.pssm
         kinase     -5P     -5G     -5A  ...      4E      4s      4t       4y
    0      AAK1  1.2242  0.4165  0.4825  ...  0.7777  0.4274  0.4274   0.4597
    1    ACVR2A  0.7057  0.8178  0.9920  ...  1.0888  1.1958  1.1958   1.0015
    ..      ...     ...     ...     ...  ...     ...     ...     ...      ...
    301    YSK4  1.0079  1.2379  1.2652  ...  0.8221  1.0770  1.0770   0.6618
    302     ZAK  1.0272  1.0890  1.1211  ...  0.6293  0.6632  0.6632   0.6943
    [303 rows x 208 columns]

    Scoring matrix
    >>> kinex.scoring_matrix
                         AAK1    ACVR2A    ACVR2B  ...      YSK1      YSK4       ZAK
    VDDEKGDSNDDYDSA -7.652375 -0.556483  0.382342  ... -4.315913 -1.380515 -2.494291
    YDSAGLLSDEDCMSV -3.490767 -0.142416  0.363835  ... -6.744147 -3.372538 -5.395429
    ..                    ...       ...       ...  ...       ...       ...       ...
    SEEEASSTEKPTKAL -3.756161  1.451276  1.884506  ... -2.819413 -0.338382 -1.325926
    ASSTEKPTKALPRKS -1.541950 -2.177326 -2.012570  ... -2.413716 -3.003581 -2.549749
    [82755 rows x 303 columns]



    """

    def __init__(self, scoring_matrix: pd.DataFrame, pssm: pd.DataFrame = get_pssm()) -> None:
        """
        Initializes the instance of the Kinex class.

        Parameters
        ----------
        pssm : pandas.DataFrame
            Normalised and scaled densiometries from PSPA experiments. 
            The table cotains on rows the kinases and on columns the positions for each aminoacid.
        scoring_matrix : pandas.DataFrame
            Table containing 82,755 experimentally identified Ser/Thr phosphosites that have been scored by 303 Ser or Thr kinase PSSM.
            The table allows the ranking of kinases, as well as the calculation of promiscuity index and median percentile for each input sequence.
        """
        self.pssm = pssm

        self.scoring_matrix = {}
        for col in scoring_matrix:
            self.scoring_matrix[col] = scoring_matrix[col].to_list()

    def __repr__(self):
        return ""

    def get_score(self, sequence: str, phospho_priming: bool = False, favorability: bool = False, method: str = 'avg') -> Score:
        """
        Checks the sequence format and it's validity.
        Computes the scores, the logarithmised scores, and ranks the kinases based on the 82,755 Ser/Thr phosphosites matrix.

        Parameters
        ----------
        sequence: str
            A string representing a peptide sequence of the 20 natural aminoacids, plus pT and pY at 9 positions surrounding a Ser/Thr phospho-acceptor.
        phospho_priming: bool, default False
            Enable/Disable the phospo-priming option. 
        favorability: bool, default False
            Considers the favorability towards either Ser or Thr phospho-acceptor.
        method: str, default 'avg'
            Accounts for multiple phosphorylation sites present in a sequence and considers either the max/min/average of the site scores. 

        Returns
        -------
        Score
            Scoring result in form of an instance of the Score class containing sequence string, kinase ranking table, median percentile and methods 
            like promiscuity_index(limit) (returns a promiscuity index) and top(number) (returns a top {number} in kinase ranking).

        Examples
        --------
        Reading necessary tables
        >>> pssm_table = pd.read_csv('./data/pssm_table.csv')
        >>> scoring_matrix = pd.read_csv('./data/scoring_matrix.csv', index_col=0)

        Initializing Kinex class instance
        >>> kinex = Kinex(pssm=pssm_table, scoring_matrix=scoring_matrix)

        Calling get_score method and saving it in result variable
        >>> result = kinex.get_score(sequence='GRNSLs*PVQA', phospho_priming=False, favorability=False)
        >>> result
        Scoring results for ['GRNSLs*PVQA']

        Getting the sequence from scoring results
        >>> result.sequence
        ['GRNSLs*PVQA']

        Kinase ranking
        >>> result.ranking
                    score  log_score  percentile_score
        kinase                                        
        NLK      8.432319   3.075929            93.394
        JNK3    19.764808   4.304862            92.232
        CDK4    12.719801   3.669004            91.457
        SMG1     1.286346   0.363278            91.203
        JNK1    19.502893   4.285616            90.773
        ...           ...        ...               ...
        TLK2     0.056653  -4.141707            10.425
        RAF1     0.132770  -2.913003             6.101
        PASK     0.118487  -3.077197             5.953
        CDC7     0.160988  -2.634972             5.613
        TLK1     0.024576  -5.346618             4.937

        [303 rows x 3 columns]

        Median percentile
        >>> result.median_percentile
        43.292

        Promiscuity index
        >>> result.promiscuity_index(90)
        7

        Top 15 kinases
        >>> result.top(15)
                    score  log_score  percentile_score
        kinase                                        
        NLK      8.432319   3.075929            93.394
        JNK3    19.764808   4.304862            92.232
        CDK4    12.719801   3.669004            91.457
        SMG1     1.286346   0.363278            91.203
        JNK1    19.502893   4.285616            90.773
        JNK2    22.196366   4.472252            90.715
        VRK2     1.523400   0.607294            90.549
        DYRK1B  12.836596   3.682191            89.887
        DYRK2   12.663987   3.662660            89.157
        P38G    20.666865   4.369248            88.478
        P38D    14.826372   3.890094            87.288
        DYRK4   10.812945   3.434688            86.502
        CDK9     7.667641   2.938783            86.102
        DYRK3    4.571399   2.192636            86.052
        CDK7     9.678460   3.274777            86.048
        """

        # start = time.perf_counter()

        sequence_format = get_sequence_format(sequence)

        # Check method and format
        if not method in ['min', 'max', 'avg', 'all']:
            raise ValueError(
                f"Method {method} is not supported. Supported methods: 'min', 'max', 'avg', 'all'")

        if sequence_format == 'unsupported':
            raise ValueError(
                f"Sequence format is not supported. Supported formats: '*' and 'central'")

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
                    sequence[id] = sequence[id][:-1] + \
                        sequence[id][-1:].lower()

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
        df_all = []
        # Iterate through every sequence for asterisk format
        if sequence_format == "*":
            number_of_seq = len(sequence)
            for id in range(number_of_seq):
                if not phospho_priming:
                    sequence[id] = sequence[id].upper()
                if method == 'all':
                    df_all.append(score(
                        sequence=sequence[id], sequence_format=sequence_format, pssm=self.pssm, favorability=favorability))
                else:
                    if id == 0:
                        df = score(sequence=sequence[id], sequence_format=sequence_format, pssm=self.pssm, favorability=favorability) / number_of_seq if method == 'avg' else score(
                            sequence=sequence[id], sequence_format=sequence_format, pssm=self.pssm, favorability=favorability)
                    else:
                        if method == 'avg':
                            df = df.add(score(
                                sequence=sequence[id], sequence_format=sequence_format, pssm=self.pssm, favorability=favorability) / number_of_seq)
                        else:
                            df = pd.concat([df, score(
                                sequence=sequence[id], sequence_format=sequence_format, pssm=self.pssm, favorability=favorability)])
                            if method == 'min':
                                df = df.groupby(df.index).min()
                            elif method == 'max':
                                df = df.groupby(df.index).max()

        # Central format
        elif sequence_format == 'central':
            if not check_sequence(sequence, sequence_format):
                raise ValueError("Invalid sequence")
            if not phospho_priming:
                sequence = sequence.upper()
            df = score(sequence=sequence, sequence_format=sequence_format,
                       pssm=self.pssm, favorability=favorability)
            sequence = [sequence]

        if not method == 'all':
            df_all.append(df)

        for table in df_all:
            table.insert(1, "log_score", np.log2(table["score"]))

            # Compute percentiles
            percentiles = []
            for kinase in table.index:
                # Get the position of the kinase score within the 82755 scored reference sequences and divide it by 82755 to get the percentile score
                percentile = (bisect.bisect_left(
                    self.scoring_matrix[kinase], table.log_score[kinase]) + 1) * 100 / 82755
                percentiles.append(percentile)

            table.insert(2, "percentile_score", percentiles)
            table.sort_values("percentile_score",
                              ascending=False, inplace=True)

        # Debugging information
        # end = time.perf_counter()
        # logging.debug(f'{sequence}, {end-start}')
        return Score(sequence, df_all)

    def get_enrichment(self, input_sites: pd.DataFrame, fc_threshold: float = 1.5, phospho_priming: bool = False, favorability: bool = False, method: str = 'avg'):
        """
        Counts the number of up/down/unregulated phosphosite sequences. 
        Using one-sided Fisher exact test determines the kinase enrichment.
        Vulcano plot of the side displaying significant enrichment for each kinase vs the corresponding p-value. 

        Parameters
        ----------
        input_sites: pd.DataFrame
            A DataFrame containing the phosphosite sequences in the first column and the logarithmised Fold Change in the second column. 
        fc_threshold: float, default 1.5
            A given threshold to check the regulation. By default the threshold is 1.5, but it can be adjusted based on the displayed plots. 

        Returns
        -------
        Enrichment
            Enrichment analysis results in a form of an instance of the Enrichment class containing the enrichment table, table with input sites and regulation,
            total numbers of upregulated/downregulated/unregulated sites, and list of failed sites. There is also a method plot(use_adjusted_pval: bool = False)

        Examples
        --------
        Reading necessary tables
        >>> pssm_table = pd.read_csv('./data/pssm_table.csv')
        >>> scoring_matrix = pd.read_csv('./data/scoring_matrix.csv', index_col=0)
        >>> input_sites = pd.read_csv('data/fam20c_cut.csv', sep='\t')

        Initializing Kinex class instance
        >>> kinex = Kinex(pssm=pssm_table, scoring_matrix=scoring_matrix)

        Calling get_score method and saving it inside of a variable
        >>> result = kinex.get_enrichment(input_sites=input_sites, fc_threshold=1.0)
        >>> result
        Total number of upregulated sites is: 6
        Total number of downregulated sites is: 77
        Total number of unregulated sites is: 24

        Accessing the values separately
        >>> result.total_upregulated
        6
        >>> result.total_downregulated
        77
        >>> result.total_unregulated
        24

        Enrichment table
        >>> result.enrichment_table
                upregulated downregulated  unregulated upregulated_enrichment_value  ...
        kinase                                                                       
        AAK1             0           8.0          3.0                     0.472527   ...
        ACVR2A           0          15.0          6.0                     0.218935   ...
        ACVR2B           0          13.0          5.0                     0.272727   ...
        AKT1           1.0             0          1.0                    13.363636   ...
        AKT2           1.0             0          1.0                          4.6   ...
        ...            ...           ...          ...                          ...   ...
        YANK2            0          10.0          1.0                     1.205128   ...
        YANK3          1.0           3.0          1.0                          4.6   ...
        YSK1             0           1.0          1.0                     3.769231   ...
        YSK4             0           2.0          1.0                     3.769231   ...
        ZAK              0           1.0          1.0                     3.769231   ...

        [284 rows x 19 columns]


        Input sites and regulation
        >>> result.input_sites
                phosphosite   log2_fc     regulation
        0     LQVKIPSKEEEAD -0.476009    unregulated
        1    LQVKIPSKEEEsAD -0.476009         failed
        2     EGRNSLSPVQATQ  0.066476    unregulated
        3     EEEADMSSPTQRT  0.649534    unregulated
        4     RRHRNLSSTTDDE  1.099660    upregulated
        ..              ...       ...            ...
        103   QTPKDGSNKSGAE -4.634178  downregulated
        104   TESGEETDLISPP -3.158583  downregulated
        105   SHPEPQTPKDSPS  0.557604    unregulated
        106   GKLCAHSQQRQYR -2.706312  downregulated
        107   KEKVHLSDSERKM -1.168763  downregulated

        [108 rows x 3 columns]

        The list of the phosphosites that did't pass the valudation
        >>> result.failed_sites
        ['LQVKIPSKEEEsAD']
        """

        if not method in ['min', 'max', 'avg']:
            raise ValueError(
                f"Method {method} is not supported. Supported methods: 'min', 'max', 'avg'")

        # start = time.perf_counter()
        df = input_sites.copy()

        df.iloc[:,0] = df.iloc[:,0].astype(str).str.replace('(ub)', '', regex=False).str.replace(
            '(ox)', '', regex=False).str.replace('(ac)', '', regex=False).str.replace('(de)', '', regex=False)

        # Empty DataFrame to store the output
        enrichment_table = pd.DataFrame(
            columns=['kinase', 'upregulated', 'downregulated', 'unregulated'])

#         logging.debug(enrichment_table)

        total_upregulated = total_downregulated = total_unregulated = 0
        regulation_list = []
        failed_sites = []
        top15_kinases_list = []

        for id in range(len(df)):
            # Get top 15 kinases, check if site is valid
            # logging.debug(f"Scoring {df.iloc[id, 0]} : {id}/{len(df) - 1}")
            try:
                top15_kinases = self.get_score(sequence=str(
                    df.iloc[id, 0]), phospho_priming=phospho_priming, favorability=favorability, method=method).top(15).index
            except ValueError:
                # logging.warning(f"Scoring of {df.iloc[id, 0]} failed")
                failed_sites.append(df.iloc[id, 0])
                regulation_list.append('failed')
                top15_kinases_list.append("")
                continue

            regulation = ""
            if float(str(df.iloc[id, 1])) >= fc_threshold:
                regulation = "upregulated"
                total_upregulated += 1
            elif float(str(df.iloc[id, 1])) <= -fc_threshold:
                regulation = "downregulated"
                total_downregulated += 1
            elif float(str(df.iloc[id, 1])) < fc_threshold:
                regulation = "unregulated"
                total_unregulated += 1

            regulation_list.append(regulation)
            top15_kinases_list.append(",".join(top15_kinases))

            enrichment_table = pd.concat([enrichment_table, pd.DataFrame(
                {"kinase": top15_kinases, regulation: np.ones(len(top15_kinases))})]).groupby('kinase').sum(numeric_only=False).reset_index()

        # Add regulation column to input_sites table
        df.insert(2, 'regulation', regulation_list)
        df.insert(3, 'top15_kinases', top15_kinases_list)

        # TODO think about background adjustment
        # Background adjustment
        if total_unregulated == 0:
            total_unregulated = np.min(
                [total_upregulated, total_downregulated])/2

        # end = time.perf_counter()
        # logging.debug(f'{end-start}')
        # logging.debug(enrichment_table)

        return Enrichment(enrichment_table, df, failed_sites, total_upregulated, total_downregulated, total_unregulated, set(self.pssm.index))
