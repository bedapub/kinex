from setuptools import setup

setup(
    name="kinex",
    version="0.0.1",
    description="A python package to compute kinase scoring and enrichment",
    py_modules=["pssm"],
    package_dir={"": "src"},
)