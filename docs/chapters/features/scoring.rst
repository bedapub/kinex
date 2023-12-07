Get scores for a phospho-sequence
=================================

1. Score a given phospho-sequence and store the results

.. code:: python

    >>> sequence = "FVKQKASQSPQKQ"
    >>> res = kinex.get_score(sequence)
    score.Score

.. note:: 

    The supported format for a phospho-sequence\
        - ``asterisk``: F_KQKAS*QSPQK
        - ``(ph)``: FKQKAS(ph)QSPQK

2. Rank the kinases for a given phospho-sequence based on their percentile score

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

4. Get median percentile (median score of kinases for a given phospho-sequence)

.. code:: python

    >>> res.median_percentile
    63.493
    numpy.float64

5. Get promiscuity index (the number of kinases scoring above the ``90th`` percentile)

.. code:: python

    >>> res.promiscuity_index()
    30
    numpy.int64