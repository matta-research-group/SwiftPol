Getting Started
===============

Getting Started with Swiftpol
=============================

Welcome to the Swiftpol documentation! This guide will help you get started with using Swiftpol for polymer molecular dynamics simulations.

Installation
------------

To install Swiftpol, you can use `pip`. Run the following command in your terminal:

```bash
pip install swiftpol
```

Alternatively, you can clone the repository and install it locally:

```bash
git clone https://github.com/<username>/swiftpol.git
cd swiftpol
pip install .
```

Usage
-----
Build a short PEG chain using Swiftpol:

```python
from swiftpol import build
from rdkit.Chem import AllChem

build.build_polymer('AAAA',
                    ['IOCCOI'], # Iodine is used to mark reaction sites
                    AllChem.ReactionFromSmarts('[C:1]-[O:2]-[I:3].[C:4]-[O:5]-[I:6]>>[C:1]-[O:2]-[C:4].[I:3][O:5][I:6]')
)
```


