import unittest
from rdkit import Chem
from swiftpol.parameterize import charge_polymer, forcefield_with_charge_handler, forcefield_with_residualmonomer
from openff.toolkit.typing.engines.smirnoff import ForceField
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
        # Test espaloma charge scheme
        charges = charge_polymer(polymer, 'espaloma')
        self.assertIsNotNone(charges)
        # Test NAGL charge scheme
        charges = charge_polymer(polymer, 'NAGL')
        self.assertIsNotNone(charges)
        # Test invalid charge scheme
        with self.assertRaises(AttributeError):
            charges = charge_polymer(polymer, 'invalid')
        
        #Test with polymer build using build
        sequence = 'ABAB'
        monomer_a_smiles = 'OC(=O)CO'
        monomer_b_smiles = 'C[C@@H](C(=O)[OH])O'
        reaction = AllChem.ReactionFromSmarts('[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]')
        polymer, A_ratio, B_ratio = build.build_linear_copolymer(sequence, monomer_a_smiles, monomer_b_smiles, reaction)
        # Test AM1_BCC charge scheme
        charges = charge_polymer(polymer, 'AM1_BCC')
        self.assertIsNotNone(charges)
        # Test espaloma charge scheme
        charges = charge_polymer(polymer, 'espaloma')
        self.assertIsNotNone(charges)
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

    def test_forcefield_with_residualmonomer(self):
        # Create a test system
        system = type('', (), {})()  # create a simple object
        system.residual_monomer = 1
        # Test with residual monomer
        forcefield = forcefield_with_residualmonomer(system, 'NAGL', ForceField("openff-2.2.0.offxml"))
        self.assertIsNotNone(forcefield)
        # Test without residual monomer
        system.residual_monomer = 0
        with self.assertRaises(AttributeError):
            forcefield = forcefield_with_residualmonomer(system, 'NAGL', ForceField("openff-2.2.0.offxml"))

if __name__ == '__main__':
    unittest.main()