import unittest
from rdkit import Chem
from swiftpol.parameterize import charge_polymer, forcefield_with_charge_handler, charge_openff_polymer
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.topology import Molecule
from swiftpol import build 
from swiftpol import demo
from rdkit.Chem import AllChem

class TestParameterize(unittest.TestCase):
    def test_charge_polymer(self):
        # Create a test polymer
        polymer = Chem.MolFromSmiles('C[C@@H](C(=O)[OH])O')
        # Test AM1_BCC charge scheme
        charges = charge_polymer(polymer, 'AM1_BCC')
        self.assertIsNotNone(charges)
        # Test espaloma charge scheme if espaloma is installed
        try:
            charges = charge_polymer(polymer, 'espaloma')
            self.assertIsNotNone(charges)
        except ImportError:
            pass  # Skip the test if espaloma is not installed
        # Test NAGL charge scheme
        charges = charge_polymer(polymer, 'NAGL')
        self.assertIsNotNone(charges)
        # Test invalid charge scheme
        with self.assertRaises(AttributeError):
            charges = charge_polymer(polymer, 'invalid')

        # Create a test polymer
        openff_polymer = Molecule.from_rdkit(polymer)
        # Test AM1_BCC charge scheme
        charges = charge_openff_polymer(openff_polymer, 'AM1_BCC')
        self.assertIsNotNone(charges)
        # Test espaloma charge scheme if espaloma is installed
        try:
            charges = charge_openff_polymer(openff_polymer, 'espaloma')
            self.assertIsNotNone(charges)
        except ImportError:
            pass  # Skip the test if espaloma is not installed
        # Test NAGL charge scheme
        charges = charge_openff_polymer(openff_polymer, 'NAGL')
        self.assertIsNotNone(charges)
        # Test invalid charge scheme
        with self.assertRaises(AttributeError):
            charges = charge_openff_polymer(openff_polymer, 'invalid')        
        
        #Test with polymer build using build
        sequence = 'ABAB'
        monomer_a_smiles = 'OC(=O)CO'
        monomer_b_smiles = 'C[C@@H](C(=O)[OH])O'
        reaction = AllChem.ReactionFromSmarts('[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]')
        polymer, A_ratio, B_ratio = build.build_linear_copolymer(sequence, monomer_a_smiles, monomer_b_smiles, reaction)
        # Test AM1_BCC charge scheme
        charges = charge_polymer(polymer, 'AM1_BCC')
        self.assertIsNotNone(charges)
        # Test espaloma charge scheme if espaloma is installed
        try:
            import espaloma
            charges = charge_polymer(polymer, 'espaloma')
            self.assertIsNotNone(charges)
        except ImportError:
            pass  # Skip the test if espaloma is not installed
        # Test NAGL charge scheme
        charges = charge_polymer(polymer, 'NAGL')
        self.assertIsNotNone(charges)

    def test_forcefield_with_charge_handler(self):
        # Create a test molecule
        molecule = Chem.MolFromSmiles('C[C@@H](C(=O)[OH])O')
        # Test with ensemble=False
        forcefield = forcefield_with_charge_handler(molecule, 'NAGL', ensemble=False)
        self.assertIsNotNone(forcefield)
        # Test with ensemble=True
        ensemble_object = demo.PLGA_system(80, 50, 1.0, 'ester', 5)
        forcefield = forcefield_with_charge_handler(ensemble_object, 'NAGL', ensemble=True)
        self.assertIsNotNone(forcefield)
        # Test with ensemble likely to contain duplicates
        ensemble_object = demo.PLGA_system(80, 50, 1.0, 'ester', 100)
        forcefield = forcefield_with_charge_handler(ensemble_object, 'NAGL', ensemble=True)
        self.assertIsNotNone(forcefield)

if __name__ == '__main__':
    unittest.main()