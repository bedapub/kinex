import bisect
from collections import namedtuple
from functools import reduce
import tomli

import numpy as np
import pandas as pd

from kinex.enrichment import Enrichment
from kinex.functions import download_file_to_resource
from kinex.resources import (
    get_pssm_ser_thr,
    get_pssm_tyr,
    get_scoring_matrix_ser_thr,
    get_scoring_matrix_tyr,
)
from kinex.score import Score
from kinex.enrichment import Enrichment
from kinex.sequence import get_sequence_object, SequenceType


EnrichmentResults = namedtuple("EnrichmentResults", ["ser_thr", "tyr", "failed_sites"])


# Load the pyproject.toml file
with open("pyproject.toml", "rb") as f:
    config = tomli.load(f)


class Kinex:
    """
    The class representing a PSSM table and a scoring matrix needed for scoring and enrichment analysis.

    Attributes
    ----------
    pssm_ser_thr : pandas.DataFrame
        Normalised and scaled densiometries from PSPA experiments.
        The table cotains on rows the kinases and on columns the positions for each aminoacid.
    scoring_matrix_ser_thr : pandas.DataFrame
        Table containing 82,755 experimentally identified Ser/Thr phosphosites that have been scored by 303 Ser or Thr kinase PSSM.
        The table allows the ranking of kinases, as well as the calculation of promiscuity index and median percentile for each input validation.

    Methods
    -------
    get_score(self, validation: str, phospho_priming: bool = False, favorability: bool = False, method: str = 'avg') -> score.Score
        Checks the validation format and it's validity; computes the scores and ranks the kinases.
    get_enrichment(self, input_sites: pd.DataFrame, fc_threshold: float = 1.5) -> enrichment.Enrichment
        Checks regulation and plots the the enrichment vs p-value.

    Examples
    --------
    PSSM table
    >>> kinex.pssm_ser_thr
         kinase     -5P     -5G     -5A  ...      4E      4s      4t       4y
    0      AAK1  1.2242  0.4165  0.4825  ...  0.7777  0.4274  0.4274   0.4597
    1    ACVR2A  0.7057  0.8178  0.9920  ...  1.0888  1.1958  1.1958   1.0015
    ..      ...     ...     ...     ...  ...     ...     ...     ...      ...
    301    YSK4  1.0079  1.2379  1.2652  ...  0.8221  1.0770  1.0770   0.6618
    302     ZAK  1.0272  1.0890  1.1211  ...  0.6293  0.6632  0.6632   0.6943
    [303 rows x 208 columns]

    Scoring matrix
    >>> kinex.scoring_matrix_ser_thr
                         AAK1    ACVR2A    ACVR2B  ...      YSK1      YSK4       ZAK
    VDDEKGDSNDDYDSA -7.652375 -0.556483  0.382342  ... -4.315913 -1.380515 -2.494291
    YDSAGLLSDEDCMSV -3.490767 -0.142416  0.363835  ... -6.744147 -3.372538 -5.395429
    ..                    ...       ...       ...  ...       ...       ...       ...
    SEEEASSTEKPTKAL -3.756161  1.451276  1.884506  ... -2.819413 -0.338382 -1.325926
    ASSTEKPTKALPRKS -1.541950 -2.177326 -2.012570  ... -2.413716 -3.003581 -2.549749
    [82755 rows x 303 columns]



    """

    def __init__(
        self,
        scoring_matrix_ser_thr: pd.DataFrame = None,
        scoring_matrix_tyr: pd.DataFrame = None,
        pssm_ser_thr: pd.DataFrame = get_pssm_ser_thr(),
        pssm_tyr: pd.DataFrame = get_pssm_tyr(),
    ) -> None:
        """
        Initializes the instance of the Kinex class.

        Parameters
        ----------
        pssm_ser_thr : pandas.DataFrame
            Normalized and scaled densiometries from PSPA experiments.
            The table cotains on rows the kinases and on columns the positions for each aminoacid.
        scoring_matrix_ser_thr : pandas.DataFrame
            Table containing 82,755 experimentally identified Ser/Thr phosphosites that have been scored by 303 Ser or Thr kinase PSSM.
            The table allows the ranking of kinases, as well as the calculation of promiscuity index and median percentile for each input validation.
        """

        # Matrix is not provided
        if scoring_matrix_ser_thr is None:
            # Trying to look for the matrix in the resources
            scoring_matrix_ser_thr = get_scoring_matrix_ser_thr()
            # Matrix is not provided and not found in the resources, download the default matrix
            if scoring_matrix_ser_thr is None:
                scoring_matrix_ser_thr_url = config["project"]["urls"][
                    "scoring_matrix_ser_thr"
                ]
                download_file_to_resource(
                    scoring_matrix_ser_thr_url, "default_scoring_matrix_ser_thr.csv.gz"
                )
                scoring_matrix_ser_thr = get_scoring_matrix_ser_thr()

        if scoring_matrix_tyr is None:
            scoring_matrix_tyr = get_scoring_matrix_tyr()
            if scoring_matrix_tyr is None:
                scoring_matrix_tyr_url = config["project"]["urls"]["scoring_matrix_tyr"]
                download_file_to_resource(
                    scoring_matrix_tyr_url, "default_scoring_matrix_tyr.csv.gz"
                )
                scoring_matrix_tyr = get_scoring_matrix_tyr()

        self.pssm_ser_thr = pssm_ser_thr
        self.pssm_tyr = pssm_tyr

        self.scoring_matrix_ser_thr = {
            col: scoring_matrix_ser_thr[col].to_list() for col in scoring_matrix_ser_thr
        }
        self.scoring_matrix_ser_thr_length = len(scoring_matrix_ser_thr)

        self.scoring_matrix_tyr = {
            col: scoring_matrix_tyr[col].to_list() for col in scoring_matrix_tyr
        }
        self.scoring_matrix_tyr_length = len(scoring_matrix_tyr)

    def __repr__(self):
        return ""

    def get_score(
        self,
        sequence: str,
        phospho_priming: bool = False,
        favorability: bool = False,
        method: str = "avg",
    ) -> Score:
        """
        Checks the validation format and it's validity.
        Computes the scores, the logarithmic scores, and ranks the kinases based on the 82,755 Ser/Thr phosphosites matrix.

        Parameters
        ----------
        sequence: str
            A string representing a peptide validation of the 20 natural aminoacids, plus pT and pY at 9 positions surrounding a Ser/Thr phospho-acceptor.
        phospho_priming: bool, default False
            Enable/Disable the phospo-priming option.
        favorability: bool, default False
            Considers the favorability towards either Ser or Thr phospho-acceptor.
        method: str, default 'avg'
            Accounts for multiple phosphorylation sites present in a validation and considers either the max/min/average of the site scores.

        Returns
        -------
        Score
            Scoring result in form of an instance of the Score class containing validation string, kinase ranking table, median percentile and methods
            like promiscuity_index(limit) (returns a promiscuity index) and top(number) (returns a top {number} in kinase ranking).

        Examples
        --------
        Reading necessary tables
        >>> pssm_table = pd.read_csv('resources/pssm_table_ser_thr.csv')
        >>> scoring_matrix = pd.read_csv('kinex/resources/scoring_matrix.csv', index_col=0)

        Initializing Kinex class instance
        >>> kinex = Kinex(pssm_ser_thr=pssm_table, scoring_matrix_ser_thr=scoring_matrix)

        Calling get_score method and saving it in result variable
        >>> result = kinex.get_score(validation='GRNSLs*PVQA', phospho_priming=False, favorability=False)
        >>> result
        Scoring results for ['GRNSLs*PVQA']

        Getting the validation from scoring results
        >>> result.validation
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

        if len(sequence) < 3:
            raise ValueError("Invalid sequence")

        if method not in ["min", "max", "avg", "all"]:
            raise ValueError(
                f"Method {method} is not supported. Supported methods: 'min', 'max', 'avg', 'all'"
            )

        if not phospho_priming:
            sequence = sequence.upper()
            
        sequence_object = get_sequence_object(sequence)
        sequence_object.preprocess_sequence()
        sequence_object.validate_sequence()

        match sequence_object.sequence_type:
            case SequenceType.SER_THR:
                pssm = self.pssm_ser_thr
                scoring_matrix = self.scoring_matrix_ser_thr
                scoring_matrix_len = self.scoring_matrix_ser_thr_length
            case SequenceType.TYR:
                pssm = self.pssm_tyr
                scoring_matrix = self.scoring_matrix_tyr
                scoring_matrix_len = self.scoring_matrix_tyr_length
            case _:
                raise ValueError("Invalid sequence type")

        scores_results = sequence_object.get_sequence_scores(pssm, favorability)

        match method:
            case "min":
                scores_results = pd.concat(scores_results)
                scores_results = [scores_results.groupby(scores_results.index).min()]
            case "max":
                scores_results = pd.concat(scores_results)
                scores_results = [scores_results.groupby(scores_results.index).max()]
            case "avg":
                scores_results = [
                    reduce(lambda df1, df2: df1.add(df2, fill_value=0), scores_results)
                    / len(scores_results)
                ]
            case "all":
                pass

        for table in scores_results:
            table.insert(1, "log_score", np.log2(table["score"]))

            percentiles = []
            for kinase in table.index:
                percentile = (
                    bisect.bisect_left(scoring_matrix[kinase], table.log_score[kinase])
                    + 1
                ) * (100 / scoring_matrix_len)
                percentiles.append(percentile)

            table.insert(2, "percentile_score", percentiles)
            table.sort_values("percentile_score", ascending=False, inplace=True)

        return Score(sequence_object, scores_results)

    def get_enrichment(
        self,
        input_sites: pd.DataFrame,
        fc_threshold: float = 1.5,
        phospho_priming: bool = False,
        favorability: bool = False,
        method: str = "avg",
    ):
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
        >>> pssm_table = pd.read_csv('resources/pssm_table_ser_thr.csv')
        >>> scoring_matrix = pd.read_csv('kinex/resources/scoring_matrix.csv', index_col=0)
        >>> input_sites = pd.read_csv('kinex/resources/fam20c_cut.csv', sep='\t')

        Initializing Kinex class instance
        >>> kinex = Kinex(pssm_ser_thr=pssm_table, scoring_matrix_ser_thr=scoring_matrix)

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

        The list of the phosphosites that didn't pass the validation
        >>> result.failed_sites
        ['LQVKIPSKEEEsAD']
        """

        if method not in ["min", "max", "avg"]:
            raise ValueError(
                f"Method {method} is not supported. Supported methods: 'min', 'max', 'avg'"
            )

        input_sites_copy = input_sites.copy(deep=True)
        input_sites_copy.iloc[:, 0] = (
            input_sites_copy.iloc[:, 0]
            .astype(str)
            .str.replace("(ub)", "", regex=False)
            .str.replace("(ox)", "", regex=False)
            .str.replace("(ac)", "", regex=False)
            .str.replace("(de)", "", regex=False)
        )

        ser_thr_enrichment = Enrichment(
            SequenceType.SER_THR, set(self.pssm_ser_thr.index)
        )
        tyr_enrichment = Enrichment(SequenceType.TYR, set(self.pssm_tyr.index))

        failed_sites = []

        for id in range(len(input_sites_copy)):
            try:
                score_result = self.get_score(
                    sequence=str(input_sites_copy.iloc[id, 0]),
                    phospho_priming=phospho_priming,
                    favorability=favorability,
                    method=method,
                )
            except ValueError:
                failed_sites.append(input_sites_copy.iloc[id, 0])
                continue

            match score_result.sequence.sequence_type:
                case SequenceType.SER_THR:
                    enrichment_object = ser_thr_enrichment
                case SequenceType.TYR:
                    enrichment_object = tyr_enrichment
                case _:
                
                    failed_sites.append(input_sites_copy.iloc[id, 0])
                    continue

            regulation = ""
            if float(str(input_sites_copy.iloc[id, 1])) >= fc_threshold:
                regulation = "upregulated"
                enrichment_object.total_upregulated += 1
            elif float(str(input_sites_copy.iloc[id, 1])) <= -fc_threshold:
                regulation = "downregulated"
                enrichment_object.total_downregulated += 1
            elif float(str(input_sites_copy.iloc[id, 1])) < fc_threshold:
                regulation = "unregulated"
                enrichment_object.total_unregulated += 1

            top15_kinases = score_result.top(15)[0].index
            enrichment_object.regulation_list.append(regulation)
            enrichment_object.top15_kinases_list.append(",".join(top15_kinases))

            # To avoid the warning here:
            # Enrichment table sometimes is empty
            # Concatenation with empty or all-NA entries is deprecated
            enrichment_object.enrichment_table = (
                pd.concat(
                    [
                        enrichment_object.enrichment_table,
                        pd.DataFrame(
                            {
                                "kinase": top15_kinases,
                                regulation: np.ones(len(top15_kinases)),
                            }
                        ),
                    ]
                )
                .groupby("kinase")
                .sum(numeric_only=False)
                .reset_index()
            )

        # Background adjustment
        ser_thr_enrichment.adjust_background_sites()
        tyr_enrichment.adjust_background_sites()

        ser_thr_enrichment.fisher_statistics()
        tyr_enrichment.fisher_statistics()

        return EnrichmentResults(ser_thr_enrichment, tyr_enrichment, failed_sites)
