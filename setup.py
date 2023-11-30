from setuptools import setup

setup(
    name="kinex",
    version="1.0.0",
    description="A python package to compute kinase scoring and enrichment",
    py_modules=["kinex"],
    package_dir={"": "src"},
    url="https://github.com/alexandra-valeanu/kinex",
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
