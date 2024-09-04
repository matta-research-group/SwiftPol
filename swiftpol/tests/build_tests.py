import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from swiftpol import build 
from rdkit.Chem.Descriptors import ExactMolWt
from openff.toolkit.topology import Molecule
import numpy as np  # Added for numpy array type checking
from openff.toolkit.topology import Topology  # Added for topology type checking

# Define all test cases

class TestBuildPLGARing(unittest.TestCase):
    def test_build_PLGA_ring(self):
        sequence = 'LLGG'
        reaction = AllChem.ReactionFromSmarts('[I:1][O:2].[I:3][C:4]>>[C:4][O:2].[I:3][I:1]')
        terminal = 'hydroxyl'
        polymer, LA_ratio, GA_ratio = build.build_PLGA_ring(sequence, reaction, terminal)
        
        # Check if the returned polymer is a RDKit Mol object
        self.assertTrue(isinstance(polymer, Chem.rdchem.Mol))
        
        # Check if the returned ratios are as expected
        self.assertEqual(LA_ratio, 50.0)
        self.assertEqual(GA_ratio, 50.0)

class TestBuildLinearCopolymer(unittest.TestCase):
    def test_build_linear_copolymer(self):
        sequence = 'ABAB'
        monomer_a_smiles = 'CC(=O)O'
        monomer_b_smiles = 'CC(=O)OC'
        reaction = AllChem.ReactionFromSmarts('[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]')
        polymer, A_ratio, B_ratio = build.build_linear_copolymer(sequence, monomer_a_smiles, monomer_b_smiles, reaction)
        
        # Check if the returned polymer is a RDKit Mol object
        self.assertTrue(isinstance(polymer, Chem.rdchem.Mol))
        
        # Check if the returned ratios are as expected
        self.assertEqual(A_ratio, 50.0)
        self.assertEqual(B_ratio, 50.0)

class TestPDI(unittest.TestCase):
    def test_PDI(self):
        # Create some RDKit molecule objects
        mol1 = Chem.MolFromSmiles('CC(=O)O')
        mol2 = Chem.MolFromSmiles('CC(=O)OC')
        mol3 = Chem.MolFromSmiles('CC(=O)OCC')
        
        chains = [mol1, mol2, mol3]
        
        pdi, mn, mw = build.PDI(chains)
        
        # Check if the returned PDI, Mn, and Mw are floats
        self.assertTrue(isinstance(pdi, float))
        self.assertTrue(isinstance(mn, float))
        self.assertTrue(isinstance(mw, float))

        # Check if the returned PDI, Mn, and Mw are as expected
        self.assertAlmostEqual(pdi, 1.0, places=2)
        self.assertAlmostEqual(mn, 60.0, places=2)
        self.assertAlmostEqual(mw, 60.0, places=2)

class TestBlockinessCalc(unittest.TestCase):
    def test_blockiness_calc(self):
        sequence = 'GLGLGLGL'
        blockiness, block_length_G, block_length_L = build.blockiness_calc(sequence)
        
        # Check if the returned blockiness, block_length_G, and block_length_L are floats
        self.assertTrue(isinstance(blockiness, float))
        self.assertTrue(isinstance(block_length_G, float))
        self.assertTrue(isinstance(block_length_L, float))

        # Check if the returned blockiness, block_length_G, and block_length_L are as expected
        self.assertAlmostEqual(blockiness, 0.0, places=2)
        self.assertAlmostEqual(block_length_G, 1.0, places=2)
        self.assertAlmostEqual(block_length_L, 1.0, places=2)

class TestCalculateBoxComponents(unittest.TestCase):
    def test_calculate_box_components(self):
        # Create some RDKit molecule objects
        mol1 = Molecule.from_smiles('CC(=O)O')
        mol2 = Molecule.from_smiles('CC(=O)OC')
        mol3 = Molecule.from_smiles('CC(=O)OCC')
        
        chains = [mol1, mol2, mol3]
        sequence = 'GLGLGLGL'
        salt_concentration = 0.1
        residual_monomer = 0.00
        
        molecules, number_of_copies, topology, box_vectors = build.calculate_box_components(chains, sequence, salt_concentration, residual_monomer)
        
        # Check if the returned molecules is a list of Molecule objects
        self.assertTrue(all(isinstance(mol, Molecule) for mol in molecules))
        
        # Check if the returned number_of_copies is a list of integers
        self.assertTrue(all(isinstance(num, int) for num in number_of_copies))
        
        # Check if the returned topology is a Topology object
        self.assertTrue(isinstance(topology, Topology))
        
        # Check if the returned box_vectors is a numpy array
        self.assertTrue(isinstance(box_vectors, np.ndarray))

class TestPLGASystem(unittest.TestCase):
    def test_PLGA_system(self):
        # Instantiate the PLGA_system class
        plga_system = build.PLGA_system(perc_lactide_target=50, length_target=10, blockiness_target=0.5, terminals='OH', num_chains=5)
        
        # Check if the PLGA_system object has been created
        self.assertIsInstance(plga_system, build.PLGA_system)
        
        # Call the charge_system method
        plga_system.charge_system()
        
        # Check if the chains have been charged
        for chain in plga_system.chains:
            self.assertIsNotNone(chain.partial_charges)
        
        # Call the build_system method
        solvated_system = plga_system.build_system(resid_monomer=0.1, salt_concentration=0.1)
        
        # Check if the solvated_system object has been created
        self.assertIsNotNone(solvated_system)

# Run all the tests

if __name__ == '__main__':
    unittest.main()
