Kinases inference analysis
==========================

1. Read your input file

.. note:: 

    Make sure to have the phospho-sequences on the first column and the log2 transformed Fold Change on the second column.

.. code:: python

    input_sites = pd.read_csv('path/to/your/input_sites.csv')
    input_sites

.. code-block:: text

                  Sequence  Fold Change: a/a' KO Clone A vs WT
    0     KLEEKQKs*DAEEDGV                          -88.159789
    1     EEDGVTGs*QDEEDSK                          -88.159789
    ..                 ...                                 ...
    462   AKEESEEs*DEDMGFG                           19.421218
    463   RNGPRDAs*PPGSEPE                           63.187703

    [464 rows x 2 columns]
    pandas.DataFrame

.. note::

    Data: CK2 catalytic sub-units knockdown


2. Run enrichment analysis with your input phospho-sequences

.. note:: 
    
    Supported methods are ``min``, ``max``, ``avg``, ``all``

.. code:: python

    enrich = kinex.get_enrichment(input_sites, fc_threshold=1.5, phospho_priming=False, favorability=True, method="max")

3. Access the total number of up-regulated, down-regulated, and un-regulated phospho-sequences

.. code:: python

    print("Total upregulated Ser/Thr kinases:", enrich.ser_thr.total_upregulated)
    print("Total downregulated Ser/Thr kinases:", enrich.ser_thr.total_downregulated)
    print("Total unregulated Ser/Thr kinases:", enrich.ser_thr.total_unregulated)

.. code-block:: text

    Total upregulated Ser/Thr kinases: 63
    Total downregulated Ser/Thr kinases: 86
    Total unregulated Ser/Thr kinases: 309

4. Check the sites that were marked as failed

.. code:: python

    enrich.failed_sites

5. Show enrichment table

.. code:: python

    enrich.ser_thr.enrichment_table

.. code-block:: text

            upregulated  downregulated  ... dominant_enrichment_value_log2 dominant_p_value_log10_abs
    kinase                                                                      
    AAK1             0            1.0   ...                      -0.263034                   0.202666
    ACVR2A        12.0           23.0   ...                      -1.562107                   3.346702
    ...            ...            ...   ...                            ...                        ...
    YSK4             0            2.0   ...                      -1.869777                    0.68218 
    ZAK            1.0            3.0   ...                       -3.47671                     1.4713
    
    [303 rows x 19 columns]
    pandas.DataFrame

6. Vulcano plot of Enrichment Odds Ratio (EOR) and p-value

.. note::

    Kinases are represented with colours corresponding to their class. 
    
.. code:: python

    fig = enrich.ser_thr.plot(use_adjusted_pval=False)
    fig.show()


.. raw:: html
    :file: ../../figures/kinase_inference.html


.. note::

    You can update your figure (marker point, axis, legend, etc.) using Plotly’s functions:
    `https://plotly.com/python/creating-and-updating-figures <https://plotly.com/python/creating-and-updating-figures>`_


7. Save the figure in a desired format

- ``.html``

.. code:: python
    
    fig.write_html("path/to/file.html")

- ``.svg``

.. code:: python

    fig.write_image("images/fig1.svg")

- ``.pdf``

.. code:: python

    fig.write_image("images/fig1.pdf")

- ``.png``

.. code:: python

    fig.write_image("images/fig1.png", scale=10)

- ``.jpeg``

.. code:: python

    fig.write_image("images/fig1.jpeg", scale=10)
