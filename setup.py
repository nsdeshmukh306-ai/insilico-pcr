from setuptools import setup, find_packages

setup(
    name            = "insilico_pcr",
    version         = "1.0.0",
    description     = "Research-grade in-silico PCR simulation with thermodynamic scoring",
    author          = "Niraj",
    python_requires = ">=3.10",
    packages        = find_packages(),
    install_requires = [
        "biopython>=1.80",
        "numpy>=1.24",
        "scipy>=1.10",
    ],
    entry_points = {
        "console_scripts": [
            "insilico_pcr = insilico_pcr.cli:main",
        ],
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
