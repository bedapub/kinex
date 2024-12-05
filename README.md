# Kinex - Kinome Exploration Tool

**Kinex** is a Python package for inferring causal kinases from phosphoproteomics data.

Paper: Kinex infers causal kinases from phosphoproteomics data. [https://doi.org/10.1101/2023.11.23.568445](https://doi.org/10.1101/2023.11.23.568445)

## Main Features

- Substrate Sequence Scoring
- Causal Kinases Inference
- Comparison with Drug Collection

## Requirements

- [conda](https://docs.conda.io/en/latest/miniconda.html)
- Python 3.11

## Installation

### From Conda

```sh
# Create and activate your conda environment
conda create --name kinex
conda activate kinex

# Install kinex package
conda install -c bioconda kinex
```

### From Source

```sh
# Create and activate a Python 3.11 conda environment 
conda create --name kinex
conda activate kinex
conda install python=3.11

# Download the package:
git clone git@github.com:bedapub/kinex.git
cd kinex

# Install the package
pip install .
```

## Quick Start

#### 1. Import Package and Create Kinex Object

```python
from kinex import Kinex
import pandas as pd
```

##### Create Kinex Object

1. With Predefined Matrices:

    ```python
    kinex = Kinex()
    ```

2. With Your Custom Matrices:

    ```python
    kinex = Kinex(scoring_matrix_ser_thr=pd.read_csv('path_to_ser_thr_matrix.csv'), scoring_matrix_tyr=pd.read_csv('path_to_tyr_matrix.csv'))
    ```

Predefined matrices can be found here:
- [Scoring Matrix for Serine/Threonine](https://zenodo.org/records/13964893/files/scoring_matrix_ser_thr_82k_sorted.csv.gz?download=1)
- [Scoring Matrix for Tyrosine](https://zenodo.org/records/13964893/files/scoring_matrix_tyr_7k_sorted.csv.gz?download=1)

#### 2. Score a Sequence

```python
sequence = "FVKQKAY*QSPQKQ"
res = kinex.get_score(sequence)
```

#### 3. Enrichment Analysis

```python
enrich = kinex.get_enrichment(input_sites, fc_threshold=1.5, phospho_priming=False, favorability=True, method="max")

enrich.ser_thr.plot()
enrich.tyr.plot()
```

## Documentation

You can find detailed documentation describing every feature of the package with examples and tutorials [here](https://kinex.readthedocs.io/en/latest/).
