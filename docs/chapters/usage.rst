Using kinex
===========

Import package and initialise data object 
-----------------------------------------

1. Import kinex

.. code:: python

	>>> from kinex import Kinex

2. Read the scoring matrix

.. code:: python

    >>> scoring_matrix = pd.read_csv("https://zenodo.org/records/10201142/files/kinex_scoring_matrix_82k_sorted.csv.gzip?download=1", compression="gzip")
    >>> scoring_matrix

                AAK1    ACVR2A  ...      YSK4       ZAK
    0     -11.147481 -6.325340  ... -6.723077 -7.402360
    1     -10.421859 -6.178601  ... -6.343452 -7.373478
    ...          ...       ...  ...
    82753   8.074270  7.289390  ...  4.525527  4.837377
    82754   8.623180  7.871226  ...  4.869195  5.062391

    [82755 rows x 303 columns]

.. note::

    You can optionally save the scoring matrix locally for future uses.

    .. code:: python

        >>> scoring_matrix.to_csv("path/to/scoring_matrix.csv")

    Or just download using the link: `https://zenodo.org/records/10201142/files/kinex_scoring_matrix_82k_sorted.csv.gzip?download=1 <https://zenodo.org/records/10201142/files/kinex_scoring_matrix_82k_sorted.csv.gzip?download=1>`_

3.  Create a kinex object

.. code:: python

    >>> kinex = Kinex(scoring_matrix)
    >>> type(kinex)
    kinex.Kinex

Use one of three main features of kinex 
-----------------------------------------

.. toctree::
   :maxdepth: 3

   features/scoring
   features/enrichment
   features/comparison

