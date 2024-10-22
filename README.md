# Kinex - Kinome Exploration Tool

**Kinex** is a Python package for infering causal kinases from phosphoproteomics data.

Paper: Kinex infers causal kinases from phosphoproteomics data. https://doi.org/10.1101/2023.11.23.568445

## Main Features

- Substrate Sequence Scoring
- Causal Kinases Inference
- Comparison with Drug Collection

## Requirements

- [conda](https://docs.conda.io/en/latest/miniconda.html)
- python 3.11

## Installation

### From Conda

```
# Create and activate your conda environment
conda create --name kinex
conda activate kinex

# Install kinex package
conda install -c bioconda kinex
```

### From source

```
# Create and activate a python 3.11 conda environment 
conda create --name kinex
conda activate kinex
conda install python=3.11

# Download the package:
git clone git@github.com:bedapub/kinex.git
cd kinex

# Install the package
pip install .
```

## Quick start

1. Import package and create Kinex object
```
from kinex import Kinex
import pandas as pd

# Read scoring matrices from zenodo
scoring_matrix_ser_thr = pd.read_csv("https://zenodo.org/records/13964893/files/scoring_matrix_ser_thr_82k_sorted.csv.gz?download=1", compression="gzip")
scoring_matrix_tyr = pd.read_csv("https://zenodo.org/records/13964893/files/scoring_matrix_tyr_7k_sorted.csv.gz?download=1", compression="gzip")

# Create Kinex object
kinex = Kinex(scoring_matrix_ser_thr, scoring_matrix_tyr)
```
2. Score a sequence
```
sequence = "FVKQKAY*QSPQKQ"
res = kinex.get_score(sequence)
```

3. Enrichment analysis
```
enrich = kinex.get_enrichment(input_sites, fc_threshold=1.5, phospho_priming=False, favorability=True, method="max")

enrich.ser_thr.plot()
enrich.tyr.plot()
```

## Documentation

You can find detailed documentation describing every feature of the package with examples and tutorials [here](https://kinex.readthedocs.io/en/latest/).
