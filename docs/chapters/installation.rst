Installing kinex
================

Requirements
------------

* `conda <https://docs.conda.io/en/latest/miniconda.html>`__
* python 3.11

.. Installation from Pip
.. ---------------------

.. 1. Create and activate a conda venv:

.. .. code:: bash

.. 	conda create --name kinex
.. 	conda activate kinex

.. 2. Install kinex from Pypi:

.. .. code:: bash

.. 	pip install kinex

Installation from Conda
------------------------

1. Create and activate your conda environment

.. code:: bash

	conda create --name kinex
	conda activate kinex

2. Install kinex package

.. code:: bash

	conda install -c bioconda kinex


Installation from Source
------------------------

1. Create and activate a python 3.11 conda environment

.. code:: bash

	conda create --name kinex
	conda activate kinex
	conda install python=3.11

2. Download the package

.. code:: bash

	git clone git@github.com:bedapub/kinex.git
	cd kinex


3. Install the package

.. code:: bash

	pip install .