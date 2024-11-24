Import Package and Initialize Kinex 
===================================

1. **Import Kinex**

.. code:: python

    from kinex import Kinex

2. **Create a Kinex Object**

- With Predefined Matrices:

It will look in the resources to find the matrices. If it doesn't find them, it will download them and save them for future use.

.. code:: python

    kinex = Kinex()

- With Your Custom Matrices:

.. code:: python

    scoring_matrix_ser_thr = pd.read_csv("path/to/scoring_matrix_ser_thr.csv")
    scoring_matrix_tyr = pd.read_csv("path/to/scoring_matrix_tyr.csv")

    kinex = Kinex(scoring_matrix_ser_thr, scoring_matrix_tyr)

.. note::

    The matrix looks like this:

    .. code:: python

                      AAK1    ACVR2A  ...      YSK4       ZAK
        0     -11.147481 -6.325340  ... -6.723077 -7.402360
        1     -10.421859 -6.178601  ... -6.343452 -7.373478
        ...          ...       ...  ...
        82753   8.074270  7.289390  ...  4.525527  4.837377
        82754   8.623180  7.871226  ...  4.869195  5.062391

.. note::

    Predefined matrices can be found here:

    - `Scoring Matrix for Serine/Threonine <https://zenodo.org/records/13964893/files/scoring_matrix_ser_thr_82k_sorted.csv.gz?download=1>`_
    - `Scoring Matrix for Tyrosine <https://zenodo.org/records/13964893/files/scoring_matrix_tyr_7k_sorted.csv.gz?download=1>`_
