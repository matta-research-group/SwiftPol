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
class TestBuildPolymer(unittest.TestCase):
    def test_build_polymer_ROP(self):
        sequence = 'AABBAABB'
        monomer_list = ['O1C(=O)C[I+][I+]OC(=O)C1', 'C[C@@H]1[I+][I+]OC(=O)[C@H](C)OC1=O'] 
        reaction = AllChem.ReactionFromSmarts('[I:1][O:2].[I:3][C:4]>>[C:4][O:2].[I:3][I:1]')
        # Test the function
        polymer = build.build_polymer(sequence = 'AAAAAGGG', 
                                        monomer_list=['O[C@H](C)C(=O)O[I]','OCC(=O)O[I]'], 
                                        reaction = AllChem.ReactionFromSmarts('[HO:1][C:2].[O:3][C:5]=[O:6]>>[C:2][O:1][C:5]=[O:6].[O:3]'),
                                        terminal ='hydroxyl')
        self.assertIsNotNone(polymer)
        # Test with an invalid sequence
        with self.assertRaises(IndexError):
            polymer = build.build_polymer('invalid', monomer_list, reaction)
        # Test with an invalid monomer list
        with self.assertRaises(IndexError):
            polymer = build.build_polymer(sequence, ['invalid'], reaction)
        # Test with an invalid reaction
        with self.assertRaises(AttributeError):
            polymer = build.build_polymer(sequence, monomer_list, 'invalid')
        
        #Test with no terminal adjustment
        polymer = build.build_polymer(sequence = 'AAAAAGGG', 
                                        monomer_list=['O[C@H](C)C(=O)O[I]','OCC(=O)O[I]'], 
                                        reaction = AllChem.ReactionFromSmarts('[HO:1][C:2].[O:3][C:5]=[O:6]>>[C:2][O:1][C:5]=[O:6].[O:3]'))
        self.assertIsNotNone(polymer)        


class TestBuildLinearCopolymer(unittest.TestCase):
    def test_build_linear_copolymer(self):
        sequence = 'ABAB'
        monomer_a_smiles = 'OC(=O)CO'
        monomer_b_smiles = 'C[C@@H](C(=O)[OH])O'
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
        mol3 = Chem.MolFromSmiles('CC(=O)OC')
        
        chains = [mol1, mol2, mol3]
        
        pdi, mn, mw = build.PDI(chains)
        
        # Check if the returned PDI, Mn, and Mw are floats
        self.assertTrue(isinstance(pdi, float))
        self.assertTrue(isinstance(mn, float))
        self.assertTrue(isinstance(mw, float))

        # Check if the returned PDI, Mn, and Mw are as expected
        self.assertAlmostEqual(pdi, 1.77, places=2)
        self.assertAlmostEqual(mn, 69.3, places=1)
        self.assertAlmostEqual(mw, 122.6, places=1)



class TestBlockinessGen(unittest.TestCase):
    def test_blockiness_Gen(self):
        sequence = 'ABABABAB'
        blockiness, block_length_A, block_length_B = build.blockiness_gen(sequence)
        
        # Check if the returned blockiness is float
        self.assertTrue(isinstance(blockiness, float))


        # Check if the returned blockiness, block_length_G, and block_length_L are as expected
        self.assertAlmostEqual(blockiness, 0.0, places=2)
        self.assertAlmostEqual(block_length_A, 1.0, places=2)
        self.assertAlmostEqual(block_length_B, 1.0, places=2)
        
        #Test case - PLA sequence
        sequence_A = 'LLLLLLLL'
        blockiness_A, block_length_A_2, block_length_B_2 = build.blockiness_gen(sequence_A)
        self.assertTrue(isinstance(blockiness_A, str))

        
#Test calculate box components
class TestCalculateBoxComponents(unittest.TestCase):
    def test_calculate_box_components(self):
        # Create a PLGA system
        x = demo.PLGA_system(75, 30, 1.7, 'ester', 2)
        x.generate_conformers()
        # Calculate box components
        molecules, number_of_copies, topology, box_vectors, residual_monomer_actual = build.calculate_box_components(chains = x.chains,
                                                                                                                    sequence=x.sequence)        

        # Check if the returned molecules is a list
        self.assertTrue(isinstance(molecules, list))
        self.assertEqual(len(molecules), 5)
        # Check if the returned number_of_copies is an integer
        self.assertTrue(isinstance(number_of_copies, list))
        print(number_of_copies)
        # Check if the returned topology is a Topology object
        self.assertTrue(isinstance(topology, Topology))
        # Check if the returned box_vectors has the correct shape
        self.assertEqual(box_vectors.shape, (3, 3))   
        # Check if the returned residual_monomer_actual is a float
        self.assertTrue(isinstance(residual_monomer_actual, float))
        # Check if the returned residual_monomer_actual is as expected
        self.assertTrue(residual_monomer_actual==0.0)     

class TestPolymerSystem(unittest.TestCase):
    def test_init(self):

        x = build.polymer_system(monomer_list=['O[C@H](C)C(=O)O[I]'], 
                    reaction = AllChem.ReactionFromSmarts('[HO:1][C:2].[O:3][C:5]=[O:6]>>[C:2][O:1][C:5]=[O:6].[O:3]'), 
                    length_target = 10, 
                    terminals = 'hydroxyl', 
                    num_chains = 5, 
                    copolymer=False)

        self.assertTrue(len(x.chains)==5)
        self.assertTrue(9 <= round(x.max_length)<= 11)
        self.assertTrue(1.5 <= round(x.PDI) <= 5)
        self.assertTrue(x.num_chains == 5)
        
        x.generate_conformers()
        self.assertTrue(len(x.chains[0].conformers[0])==len(x.chains[0].atoms))
        x.charge_system()
        self.assertTrue(len(x.chains[0].partial_charges)==len(x.chains[0].atoms))
        from openff.units import unit



# Run

if __name__ == '__main__':
    unittest.main()
