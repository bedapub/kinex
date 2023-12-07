Drug comparison
===============

1. Import Comparison class from kinex

.. code:: python

    >>> from kinex import Comparison

2. Initialize a Comparison object

.. code:: python

    >>> comp = Comparison()


Compare multiple experiments with each other
--------------------------------------------

3. Specify the path to your enrichment tables

.. code:: python

    >>> data_path = "path/to/your/tables"


.. note:: 

    The directory should have multiple ``.csv`` files that contain enrichment analysis result tables from kinex

    .. code::

        ├── tables
            ├── table0.csv
            ├── table1.csv
            ├── table2.csv
            └── table3.csv

4. Perform the comparison

.. code:: python

    >>> fig = comp.get_comparison(data_path=data_path, method='mds')

.. note:: 

    Supported methods are ``UMAP``, ``MDS`` and ``t-SNE``.

5. Show the graph

.. code:: python

    >>> fig.show()

.. raw:: html
    :file: ../../figures/comparison_multiple_drugs.html


.. note::

    You can update your figure (marker point, axis, legend, etc.) using Plotly’s functions:
    `https://plotly.com/python/creating-and-updating-figures <https://plotly.com/python/creating-and-updating-figures>`_


6. You can optionally :ref:`save the plot in a desired format<Save the plot in a desired format>`


Compare an experiment to the existing collection of drug profiles
-----------------------------------------------------------------

3. Read the enrichment analysis result table

.. code:: python

    >>> input_data = pd.read_csv('tables/table1.csv', index_col=0)

.. note::

    The table should contain ``dominant_enrichment_value_log2`` and ``dominant_p_value_log10_abs`` columns

    .. code::

        dominant_enrichment_value_log2  dominant_p_value_log10_abs  
                             0.868162                    0.821932  
                            -0.785398                    0.707911  
                            -0.934463                    0.901927  
                            -1.369094                    0.000000  
                            -1.474303                    0.000000  
                                ...                         ...  
                            -2.914661                    2.022525  
                            -2.490535                    1.691968  
                            -2.920072                    0.000000  
                            -1.551978                    0.795959  
                            -2.986266                    1.521982  

        [303 rows x 4 columns]

4. Perform the comparison

.. code:: python

    >>> fig = comp.get_comparison(input_data=input_data, method='tsne')

.. note:: 

    Supported methods are ``UMAP``, ``MDS``, and ``t-SNE``


5. Show the graph

.. note::

    Each point represents a sample, which in this context means a unique combination of drug, 
    concentration, the duration of the treatment, the cell line used, and the running index of replicate. 
    The origin point (0, 0) represents the effect of vehicle control, i.e. no changed kinase activities. 
    If you hover over each point you can see the sample's name.


.. code:: python

    >>> fig.show()

.. raw:: html
    :file: ../../figures/comparison_input.html


Save the plot in a desired format
---------------------------------

- ``.html``

.. code:: python
    
    >>> fig.write_html("path/to/file.html")

- ``.svg``

.. code:: python

    >>> fig.write_image("images/fig1.svg")

- ``.pdf``

.. code:: python

    >>> fig.write_image("images/fig1.pdf")

- ``.png``

.. code:: python

    >>> fig.write_image("images/fig1.png")

- ``.jpeg``

.. code:: python

    >>> fig.write_image("images/fig1.jpeg")
