Get scores for a phospho-sequence
=================================

1. Score a given phospho-sequence and store the results

.. code:: python

    sequence = "FVKQKASQSPQKQ"
    res = kinex.get_score(sequence)

.. note:: 

    The supported format for a phospho-sequence\
        - ``asterisk``: F_KQKAS*QSPQK
        - ``(ph)``: FKQKAS(ph)QSPQK
        - ``central``: FVKQKASQSPQKQ

2. Rank the kinases for a given phospho-sequence based on their percentile score

.. code:: python

    res.ranking

.. code-block:: text

                score  log_score  percentile_score
    kinase                                        
    PRKD2   30.071143   4.910308            98.923
    ATR      6.798658   2.765250            98.531
    MNK1     9.017093   3.172662            97.935
    ...           ...        ...               ...
    JNK3     0.062386  -4.002623             2.043
    JNK1     0.021452  -5.542765             0.453

    [303 rows x 3 columns]
    pandas.DataFrame

3. Get the top n kinases from percentile ranking

.. code:: python

    res.top(3)

.. code-block:: text

                score  log_score  percentile_score
    kinase                                        
    PRKD2   30.071143   4.910308            98.923
    ATR      6.798658   2.765250            98.531
    MNK1     9.017093   3.172662            97.935
    pandas.DataFrame

4. Get median percentile (median score of kinases for a given phospho-sequence)

.. code:: python

    res.median_percentile

.. code-block:: text
    63.493
    numpy.float64

5. Get promiscuity index (the number of kinases scoring above the ``90th`` percentile)

.. code:: python

    res.promiscuity_index()

.. code-block:: text
    30
    numpy.int64