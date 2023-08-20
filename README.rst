kinex
=====

A python package to compute kinase scoring and enrichment.

Installation
============

Run the following to install:

.. code:: bash

	pip install -e .


Using kinex
===========

Import package and initialise data object 
-----------------------------------------

1. Import kinex

.. code:: python

	>>> from kinex import Kinex

2. Read the scoring matrix

.. code:: python

    >>> scoring_matrix = pd.read_csv('path/to/your/scoring_matrix.csv', index_col=0)
    >>> scoring_matrix

                        AAK1    ACVR2A  ...      YSK4       ZAK
    VDDEKGDSNDDYDSA -7.652375 -0.556483  ... -1.380515 -2.494291
    YDSAGLLSDEDCMSV -3.490767 -0.142416  ... -3.372538 -5.395429
    ...                   ...       ...  ...       ...       ...
    SEEEASSTEKPTKAL -3.756161  1.451276  ... -0.338382 -1.325926
    ASSTEKPTKALPRKS -1.541950 -2.177326  ... -3.003581 -2.549749

    [82755 rows x 303 columns]

3.  Create a kinex object.

.. code:: python

    >>> kinex = Kinex(scoring_matrix)
    >>> type(kinex)
    kinex.Kinex


Get scores for a sequence
-------------------------

1. Score a given sequence and store the results.

.. code:: python

    >>> sequence = "FVKQKASQSPQKQ"
    >>> res = kinex.get_score(sequence)
    score.Score

2. Rank the kinases, for a given sequence, based on their percentile score.

.. code:: python

    >>> res.ranking

                score  log_score  percentile_score
    kinase                                        
    PRKD2   30.071143   4.910308            98.923
    ATR      6.798658   2.765250            98.531
    MNK1     9.017093   3.172662            97.935
    BUB1     9.066003   3.180467            97.679
    MNK2     8.275760   3.048892            97.455
    ...           ...        ...               ...
    GRK3     0.140495  -2.831414            11.024
    MTOR     0.203042  -2.300148             9.881
    CK1D     0.593436  -0.752837             6.265
    JNK3     0.062386  -4.002623             2.043
    JNK1     0.021452  -5.542765             0.453

    [303 rows x 3 columns]
    pandas.DataFrame

3. Get the top n kinases from percentile ranking

.. code:: python

    >>> n = 5
    >>> res.top(n)

                score  log_score  percentile_score
    kinase                                        
    PRKD2   30.071143   4.910308            98.923
    ATR      6.798658   2.765250            98.531
    MNK1     9.017093   3.172662            97.935
    BUB1     9.066003   3.180467            97.679
    MNK2     8.275760   3.048892            97.455
    pandas.DataFrame

4. Get median percentile. 

.. code:: python

    >>> res.median_percentile
    63.493
    numpy.float64

5. Get promiscuity index. By default, the limit is the 90th percentile.

.. code:: python

    >>> res.promiscuity_index()
    30
    numpy.int64
 

Enrichment analysis
-------------------

1. Read your input sequences file. Make sure to have on first column the sequences and on second column the logarithmised Fold Change. 

.. code:: python

    >>> input_sites = pd.read_csv('path/to/your/input_sites.csv')
    >>> input_sites
                  site   log2_fc
    0    LQVKIPSKEEEAD -0.476009
    1    EGRNSLSPVQATQ  0.066476
    ..             ...       ...
    107  GKLCAHSQQRQYR -2.706312
    108  KEKVHLSDSERKM -1.168763

    [109 rows x 2 columns]
    pandas.DataFrame

2. Run enrichment analysis with your input phosphosite sequences.

.. code:: python

    >>> enrich = kinex.get_enrichment(input_sites, fc_threshold=1, phospho_priming=False, favorability=True, method='avg')
    >>> enrich
    Total number of upregulated sites is: 3
    Total number of downregulated sites is: 70
    Total number of unregulated sites is: 34
    enrichment.Enrichment

3. Access the total number of upregulated, downregulated and unregulated sites. 

.. code:: python

    >>> enrich.total_upregulated
    3
    int
    >>> enrich.total_downregulated
    70
    int
    >>> enrich.total_unregulated
    34
    int

4. Check the sites that were marked as failed. 

.. code:: python

    >>> enrich.failed_sites
    ['PEVVGSDSEVEG', 'EEEADMIJSSPTQRT']
    list

5. Check the regulation of each sequence. 

.. code:: python

    >>> enrich.input_sites
                  site   log2_fc     regulation
    0    LQVKIPSKEEEAD -0.476009    unregulated
    1    EGRNSLSPVQATQ  0.066476    unregulated
    ..             ...       ...            ...
    107  GKLCAHSQQRQYR -2.706312  downregulated
    108  KEKVHLSDSERKM -1.168763    unregulated

    [109 rows x 3 columns]
    pandas.DataFrame

6. Show enrichment table.

.. code:: python

    >>> enrich.enrichment_table

    kinase  upregulated  downregulated  ...  dominant_p_value_log10_abs  dominant_adjusted_p_value_log10_abs
    AAK1              0            7.0  ...                         0.0                                  0.0 
    ACVR2A            0           12.0  ...                    0.021065                                  0.0 
    ...             ...            ...  ...                         ...                                  ...
    YSK4              0            2.0  ...                         0.0                                  0.0 
    ZAK               0            1.0  ...                         0.0                                  0.0 
         
    [282 rows x 19 columns]
    pandas.DataFrame

7. Vulcano plot of enrichment vs p-value. Kinases are represented with colours corresponding to their class. 

.. code:: python

    >>> fig = enrich.plot(use_adjusted_pval=False)

.. raw:: html
    :file: ../figures/fig.html

8. Save the figure with html format.

.. code:: python

    >>> fig.write_html('path/to/your/figure.html')
