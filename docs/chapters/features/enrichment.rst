Enrichment analysis
===================

1. Read your input sequences file. Make sure to have on first column the sequences and on second column the log2 transformed Fold Change.

.. code:: python

    >>> input_sites = pd.read_csv('path/to/your/input_sites.csv')
    >>> input_sites
              sequence   log2_fc
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

5. Check the regulation of each sequence and top 15 kinases most likely to target each sequence

.. code:: python

    >>> enrich.input_sites
                  site   log2_fc     regulation top15_kinases
    0    LQVKIPSKEEEAD -0.476009    unregulated CAMK2B,CAMK2G,GSK3B,FAM20C,CAMK2A,MAPKAPK2,GRK...
    1    EGRNSLSPVQATQ  0.066476    unregulated NLK,SMG1,CDK4,JNK3,DYRK1B,JNK1,DYRK2,JNK2,P38G...
    ..             ...       ...            ... ...
    107  GKLCAHSQQRQYR -2.706312  downregulated CHAK1,IRE1,IRE2,SMG1,HRI,ATR,DNAPK,CHAK2,STK33...
    108  KEKVHLSDSERKM -1.168763    unregulated CK2A1,CK2A2,DRAK1,MEK1,GRK5,ACVR2B,SMG1,GRK6,A...

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
    :file: ../../figures/kinase_inference.html


.. note::

    Data: CK2 catalytic sub-units knockdown



8. Save the figure in a desired format.


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

    >>> fig.write_image("images/fig1.png", scale=10)

- ``.jpeg``

.. code:: python

    >>> fig.write_image("images/fig1.jpeg", scale=10)
