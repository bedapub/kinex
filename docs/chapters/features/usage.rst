Import package and initialise Kinex 
===================================

1. Import kinex

.. code:: python

	from kinex import Kinex

2. Read the scoring matrix

.. code:: python

    scoring_matrix_ser_thr = pd.read_csv("https://zenodo.org/records/13964893/files/scoring_matrix_ser_thr_82k_sorted.csv.gz?download=1", compression="gzip")
    scoring_matrix_tyr = pd.read_csv("https://zenodo.org/records/13964893/files/scoring_matrix_tyr_7k_sorted.csv.gz?download=1", compression="gzip")
    scoring_matrix_ser_thr

.. code:: python

                AAK1    ACVR2A  ...      YSK4       ZAK
    0     -11.147481 -6.325340  ... -6.723077 -7.402360
    1     -10.421859 -6.178601  ... -6.343452 -7.373478
    ...          ...       ...  ...
    82753   8.074270  7.289390  ...  4.525527  4.837377
    82754   8.623180  7.871226  ...  4.869195  5.062391

    [82755 rows x 303 columns]

.. note::

    You can optionally save the scoring matrix locally for faster use in the future.

    .. code:: python

        scoring_matrix_ser_thr.to_csv("scoring_matrix_ser_thr.csv")
        scoring_matrix_tyr.to_csv("scoring_matrix_tyr.csv")

    Or just download using the links: 
    `https://zenodo.org/records/13964893/files/scoring_matrix_ser_thr_82k_sorted.csv.gz?download=1 <https://zenodo.org/records/13964893/files/scoring_matrix_ser_thr_82k_sorted.csv.gz?download=1>`_
    `https://zenodo.org/records/13964893/files/scoring_matrix_tyr_7k_sorted.csv.gz?download=1 <https://zenodo.org/records/13964893/files/scoring_matrix_tyr_7k_sorted.csv.gz?download=1>`_

3.  Create a kinex object

.. code:: python

    kinex = Kinex(scoring_matrix_ser_thr, scoring_matrix_tyr)
