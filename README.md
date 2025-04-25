![title](Repo/logo.jpg)

[![Zenodo](https://zenodo.org/badge/798303954.svg)](https://doi.org/10.5281/zenodo.15111991)
[![PyTest](https://github.com/matta-research-group/SwiftPol/actions/workflows/SwiftPol_tests.yml/badge.svg)](https://github.com/matta-research-group/SwiftPol/actions/workflows/SwiftPol_tests.yml)
[![Documentation](https://img.shields.io/badge/Documentation-Online-brightgreen)](https://matta-research-group.github.io/SwiftPol)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://github.com/yourusername/yourrepository/blob/main/LICENSE)
[![GitHub last commit](https://img.shields.io/github/last-commit/matta-research-group/SwiftPol)](https://github.com/matta-research-group/SwiftPol/commits/main)
[![GitHub stars](https://img.shields.io/github/stars/matta-research-group/SwiftPol)](https://github.com/matta-research-group/SwiftPol/stargazers)


# SwiftPol
Tools for building polydisperse polymer systems for molecular dynamics.

This repository is currently under development. To do a development install including dependencies, clone this repository and type

`cd SwiftPol`

`conda env create --file Dev_tools/swiftpol.yml`

`conda activate swiftpol`

`pip install -e .`


### Contributing to SwiftPol

Features should be developed on branches. To create and switch to a branch, use the command

`git checkout -b < new_branch_name >`

To switch to an existing branch, use

`git checkout < new_branch_name >`

To submit your feature to be incorporated to the main branch, you should submit a Pull Request. The repository maintainers will review your pull request before accepting your changes.

### Example Notebooks
Examples of using SwiftPol code to build different polymers can be found at [Example Notebooks](Example_Notebooks/)
-  [Building a PLGA system](Example_Notebooks/PLGA_demo.ipynb)
-  [Building Chitin](Example_Notebooks/Chitin.ipynb)
-  [Constructing Reaction SMARTS](Example_Notebooks/rxn_smarts.ipynb)

### OpenEye License Guidance
[Instructions for implementing an OpenEye License (not essential but speeds up conformation determination)](https://docs.eyesopen.com/toolkits/python/quickstart-python/license.html)

Citation for OpenEye OMEGA:

Hawkins, P. C. D.; Skillman, A. G.; Warren, G. L.; Ellingson, B. A.; Stahl, M. T. Conformer Generation with OMEGA: Algorithm and Validation Using High Quality Structures from the Protein Databank and Cambridge Structural Database. J. Chem. Inf. Model. 2010, 50 (4), 572–584. https://doi.org/10.1021/ci100031x.

We acknowledge the use of a free academic license provided by OpenEye, Candence Molecular Sciences, Santa Fe, NM https://www.eyesopen.com/


### File Tree
```
.
├── CODE_OF_CONDUCT.md
├── Dev_tools
│   └── swiftpol.yml
├── Example_Notebooks
│   ├── Chitin.ipynb
│   ├── Crosslinked_PEDGA.ipynb
│   ├── PLGA_demo.ipynb
│   └── rxn_smarts.ipynb
├── LICENSE
├── pyproject.toml
├── README.md
├── Repo
│   └── logo.jpg
├── setup.py
├── swiftpol
│   ├── build.py
│   ├── demo.py
│   ├── __init__.py
│   ├── parameterize.py
│   └── _version.py
└── tests
    ├── build_tests.py
    ├── demo_tests.py
    └── parameterize_tests.py
```

### Code of Conduct
All repository contributors should follow our [code of conduct](CODE_OF_CONDUCT.md)

#### Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.

Copyright (c) 2024, Hannah Turney, Matta Research Group

