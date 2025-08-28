"""A Python package for building in silico polymer systems"""

# Init
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import random
import numpy as np
import numpy as np
import time

try:
    from openeye import oechem

    oechem_imported = True
except:
    import warnings

    warnings.warn(
        "OpenEye is not installed. You will not be able to use OpenEye Toolkits for conformer generation."
    )


import openmm
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import get_data_file_path
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, AmberToolsToolkitWrapper
from openff.units import unit
from pandas import read_csv


from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box


from functools import reduce
from statistics import mean
from rdkit.Chem.Descriptors import ExactMolWt
from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box

from swiftpol import build
from swiftpol import parameterize
from swiftpol import __version__
