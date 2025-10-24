![title](Repo/logo.jpg)

[![status](https://joss.theoj.org/papers/8ddec851cf6ac8c0ae55ed4024b75b79/status.svg)](https://joss.theoj.org/papers/8ddec851cf6ac8c0ae55ed4024b75b79)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15111991.svg)](https://doi.org/10.5281/zenodo.15111991)
[![PyTest](https://github.com/matta-research-group/SwiftPol/actions/workflows/SwiftPol_tests.yml/badge.svg)](https://github.com/matta-research-group/SwiftPol/actions/workflows/SwiftPol_tests.yml)
[![Documentation](https://img.shields.io/badge/Documentation-Online-brightgreen)](https://matta-research-group.github.io/SwiftPol)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://github.com/yourusername/yourrepository/blob/main/LICENSE)
[![GitHub last commit](https://img.shields.io/github/last-commit/matta-research-group/SwiftPol)](https://github.com/matta-research-group/SwiftPol/commits/main)
[![GitHub stars](https://img.shields.io/github/stars/matta-research-group/SwiftPol)](https://github.com/matta-research-group/SwiftPol/stargazers)



# SwiftPol
Tools for building polydisperse polymer systems for molecular dynamics.

### Installation

```bash

git clone https://github.com/matta-research-group/SwiftPol

cd SwiftPol
# install requirements into new environment
conda env create --file Dev_tools/swiftpol.yml

conda activate swiftpol

pip install -e .

```

### Contributing to SwiftPol

Features should be developed on branches. To create and switch to a branch, use the command

`git checkout -b < new_branch_name >`

To switch to an existing branch, use

`git checkout < new_branch_name >`

To submit your feature to be incorporated to the main branch, you should submit a Pull Request. The repository maintainers will review your pull request before accepting your changes.

### Example Notebooks
Examples of using SwiftPol code to build different polymers can be found at [Example Notebooks](Example_Notebooks/)
-  [Building a PLGA System](Example_Notebooks/PLGA_demo.ipynb)
-  [Building Chitin](Example_Notebooks/Chitin.ipynb)
-  [Constructing Reaction SMARTS](Example_Notebooks/rxn_smarts.ipynb)
-  [Building Crosslinked PEGDA Networks](Example_Notebooks/Crosslinking_demo.ipynb)

### Dependencies
[YAML environment file listing all required packages to use SwiftPol](https://github.com/matta-research-group/SwiftPol/blob/main/Dev_tools/swiftpol.yml)

To use SwiftPol please download the following packages:
- RDkit
- openff-interchange
- openff-toolkit
- openff-nagl
- openeye-toolkits
- openff-units==0.2.2
- dgl==2.0.0

optional dependencies:
- espaloma charge


for quick dependency and swiftpol install using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), run the following commands in bash
```bash
conda create -n swiftpol python=3.10 rdkit openff-interchange openff-toolkit openff-nagl dgl=2.0.0 openeye-toolkits nglview openff-units=0.2.2 -c conda-forge -c dglteam -c openeye
conda activate swiftpol
git clone https://github.com/matta-research-group/SwiftPol
cd SwiftPol
pip install -e .
```

### Installing With a Docker Container

To install SwiftPol as a docker container, install [docker desktop](https://docs.docker.com/desktop/) and run the following command in the terminal

```bash
docker build --no-cache --platform=linux/amd64 -f Dockerfile -t my-image -t setup/swiftpol:latest .
```

to open a jupyter notebook run

```bash
docker run --rm -it -p 8888:8888 setup/swiftpol:latest
```

### Reporting Issues and Seeking Support

If you encounter issues or erroneous behaviour when using SwiftPol, or would like support with using the software tools please [raise an issue](https://github.com/matta-research-group/SwiftPol/issues).



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

### Citing SwiftPol
If you have used SwiftPol please cite our [paper](https://joss.theoj.org/papers/10.21105/joss.08053)

### Code of Conduct
All repository contributors should follow our [code of conduct](CODE_OF_CONDUCT.md)

#### Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.

Copyright (c) 2024, Hannah Turney, Matta Research Group

