import unittest
from rdkit import Chem
from swiftpol.parameterize import (
    charge_polymer,
    charge_openff_polymer,
)
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule
from swiftpol import build
from swiftpol import demo
from rdkit.Chem import AllChem


class TestParameterize(unittest.TestCase):
    def test_charge_polymer(self):
        # Create a test polymer
        polymer = Chem.MolFromSmiles("C[C@@H](C(=O)[OH])O")
        # Test AM1_BCC charge scheme
        charges = charge_polymer(polymer, "AM1_BCC")
        self.assertIsNotNone(charges)
        # Test espaloma charge scheme if espaloma is installed
        try:
            charges = charge_polymer(polymer, "espaloma")
            self.assertIsNotNone(charges)
        except ImportError:
            pass  # Skip the test if espaloma is not installed
        # Test NAGL charge scheme
        charges = charge_polymer(polymer, "NAGL")
        self.assertIsNotNone(charges)
        # Test invalid charge scheme
        with self.assertRaises(AttributeError):
            charges = charge_polymer(polymer, "invalid")

        # Create a test polymer
        openff_polymer = Molecule.from_rdkit(polymer)
        # Test AM1_BCC charge scheme
        charges = charge_openff_polymer(openff_polymer, "AM1_BCC")
        self.assertIsNotNone(charges)
        # Test espaloma charge scheme if espaloma is installed
        try:
            charges = charge_openff_polymer(openff_polymer, "espaloma")
            self.assertIsNotNone(charges)
        except ImportError:
            pass  # Skip the test if espaloma is not installed
        # Test NAGL charge scheme
        charges = charge_openff_polymer(openff_polymer, "NAGL")
        self.assertIsNotNone(charges)
        # Test invalid charge scheme
        with self.assertRaises(AttributeError):
            charges = charge_openff_polymer(openff_polymer, "invalid")

        # Test with overwrite = False
        openff_polymer2 = Molecule.from_rdkit(polymer)
        charges = charge_openff_polymer(openff_polymer, "NAGL", overwrite=False)
        self.assertIsNotNone(charges)
        self.assertIsNone(openff_polymer2.partial_charges)

        # Test with polymer build using build
        sequence = "ABAB"
        monomer_a_smiles = "OC(=O)CO"
        monomer_b_smiles = "C[C@@H](C(=O)[OH])O"
        reaction = AllChem.ReactionFromSmarts(
            "[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]"
        )
        polymer, A_ratio, B_ratio = build.build_linear_copolymer(
            sequence, monomer_a_smiles, monomer_b_smiles, reaction
        )
        # Test AM1_BCC charge scheme
        charges = charge_polymer(polymer, "AM1_BCC")
        self.assertIsNotNone(charges)
        # Test espaloma charge scheme if espaloma is installed
        try:
            import espaloma

            charges = charge_polymer(polymer, "espaloma")
            self.assertIsNotNone(charges)
        except ImportError:
            pass  # Skip the test if espaloma is not installed
        # Test NAGL charge scheme
        charges = charge_polymer(polymer, "NAGL")
        self.assertIsNotNone(charges)

if __name__ == "__main__":
    unittest.main()
