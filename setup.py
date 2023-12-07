from setuptools import setup

setup(
    name="kinex",
    version="1.0.2",
    description="A python package to compute kinase scoring and enrichment",
    packages=["kinex", "kinex.data", "kinex.tests", "kinex.tests.data"],
    package_dir={"kinex": "src",
                 "kinex.data": "src/data",
                 "kinex.tests": "tests",
                 "kinex.tests.data": "tests/data"},
    package_data={'kinex.data': ['experiments.json', 'groups.json', 'pssm_table.csv'],
                  'kinex.tests.data': ['test_scoring_matrix.csv', 'test_pssm_table.csv', 'test_input_sites.csv']},
    include_package_data=True,
    url="https://github.com/bedapub/kinex",
    author="Alexandra Valeanu",
    author_email="valeanualexandra17@gmail.com",
    python_requires='>=3.8',
    install_requires=[
        "scipy >= 1.10.0",
        "numpy >= 1.19.5",
        "nbformat>=4.2.0",
        "pandas",
        "statsmodels",
        "plotly",
        "scikit-learn",
        "umap-learn",
        "importlib-resources"
    ],
    extras_require={
        "dev": [
            "sphinx",
            "furo",
        ],
    },
)
