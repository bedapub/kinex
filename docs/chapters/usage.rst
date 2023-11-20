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

Use one of three main features of kinex 
-----------------------------------------

.. toctree::
   :maxdepth: 3

   features/scoring
   features/enrichment
   features/comparison

