import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from swiftpol import build 
from swiftpol import demo
import numpy as np  # Added for numpy array type checking
from openff.toolkit.topology import Topology  # Added for topology type checking
import warnings
warnings.filterwarnings("ignore")

# Define all test cases

class TestBuildPLGARing(unittest.TestCase):
    def test_build_PLGA_ring(self):
        sequence = 'LLGG'
        reaction = AllChem.ReactionFromSmarts('[I:1][O:2].[I:3][C:4]>>[C:4][O:2].[I:3][I:1]')
        terminal = 'hydroxyl'
        polymer, LA_ratio, GA_ratio = demo.build_PLGA_ring(sequence, reaction, terminal)
        
        # Check if the returned polymer is a RDKit Mol object
        self.assertTrue(isinstance(polymer, Chem.rdchem.Mol))
        
        # Check if the returned ratios are as expected
        self.assertEqual(LA_ratio, 50.0)
        self.assertEqual(GA_ratio, 50.0)
        

class TestBlockinessPLGA(unittest.TestCase):
    def test_blockiness_PLGA(self):
        sequence = 'GLGLGLGL'
        blockiness, block_length_G, block_length_L = demo.blockiness_PLGA(sequence)
        
        # Check if the returned blockiness is float
        self.assertTrue(isinstance(blockiness, float))


        # Check if the returned blockiness, block_length_G, and block_length_L are as expected
        self.assertAlmostEqual(blockiness, 0.0, places=2)
        self.assertAlmostEqual(block_length_G, 1.0, places=2)
        self.assertAlmostEqual(block_length_L, 1.0, places=2)
        
        #Test case - PLA sequence
        sequence_L = 'LLLLLLLL'
        blockiness_L, block_length_G_2, block_length_L_2 = demo.blockiness_PLGA(sequence_L)
        self.assertTrue(isinstance(blockiness_L, str))
        

class TestPLGABuild(unittest.TestCase):
    def test_plga_build(self):
        x = demo.PLGA_system(80, 10, 1.0, 'ester', 5)
        
        self.assertTrue(len(x.chains)==5)
        self.assertTrue(72 <= round(x.lactide_actual) <= 88)
        self.assertTrue(9 <= round(x.max_length)<= 11)
        self.assertTrue(1.5 <= round(x.PDI) <= 5)
        print(x.PDI)
        self.assertTrue(x.mean_blockiness==1)
        
        x.generate_conformers()
        self.assertTrue(len(x.chains[0].conformers[0])==len(x.chains[0].atoms))
        x.charge_system()
        self.assertTrue(len(x.chains[0].partial_charges)==len(x.chains[0].atoms))
        from openff.units import unit
        solv_system = x.solvate_system(resid_monomer = 0.5, salt_concentration = 0.1 * unit.mole / unit.liter)
        self.assertAlmostEqual(x.residual_monomer,0.5,places=1)
        self.assertIsNotNone(solv_system)  


# Run

if __name__ == '__main__':
    unittest.main()
