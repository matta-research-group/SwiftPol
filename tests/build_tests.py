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
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')   

# Define all test cases
class TestBuildPolymer(unittest.TestCase):
    from rdkit import RDLogger 
    RDLogger.DisableLog('rdApp.*')   
    def test_build_polymer(self):
        sequence = 'AABBAABB'
        monomer_list = ['OC(=O)COI', 'C[C@@H](C(=O)[OH])OI']
        reaction = AllChem.ReactionFromSmarts('[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]')
        # Test the function
        polymer = build.build_polymer(sequence = sequence, 
                                        monomer_list=monomer_list, 
                                        reaction = reaction,
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
        polymer = build.build_polymer(sequence = 'AAAAABBBBB', 
                                        monomer_list=['O[C@H](C)C(=O)O[I]','OCC(=O)O[I]'], 
                                        reaction = reaction)
        self.assertIsNotNone(polymer)
        for atom in polymer.GetAtoms():
            assert atom.GetPDBResidueInfo() is not None



        #Test with cellulose
        monomer_list = ["C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OI)O)O)OI)O"] # Glucose with iodinated 1,4 positions
        reaction = AllChem.ReactionFromSmarts('[I:1][O:2].[I:3][O:4][C:5]>>[C:5][O:2].[I:3][I:1].[O:4]')
        sequence = 'AAA'
        polymer = build.build_polymer(sequence = sequence,
                                        monomer_list = monomer_list,
                                        reaction = reaction,
                                        terminal = 'hydroxyl')
        self.assertIsNotNone(polymer)

        #Test with chitin
        monomer_list=['C([C@@H]1[C@H]([C@@H]([C@H](C(O1)OI)N)O)OI)O']
        reaction=AllChem.ReactionFromSmarts('[I:1][O:2].[I:3][O:4][C:5]>>[C:5][O:2].[I:3][I:1].[O:4]')
        terminal='hydroxyl'
        sequence = 'AAA'
        polymer = build.build_polymer(sequence = sequence,
                                        monomer_list = monomer_list,
                                        reaction = reaction,
                                        terminal = terminal)
        self.assertIsNotNone(polymer)

        #Test with polyethylene
        monomer_list = ['ICCI']
        reaction = AllChem.ReactionFromSmarts('[C:1][C:2][I:3].[C:4][C:5][I:6]>>[C:1][C:2][C:5][C:4].[I:3][I:6]')
        terminal='hydroxyl'
        sequence = 'AAAAAAAA'
        polymer = build.build_polymer(sequence = sequence,
                                        monomer_list = monomer_list,
                                        reaction = reaction,
                                        terminal = terminal)
        self.assertIsNotNone(polymer)

        #Test with polythiophene
        monomer_list = ['S1C(CCI)=CC=C(CCI)1']
        reaction = AllChem.ReactionFromSmarts('[H:7]-[C:1]-[C:2]-[I:3].[H:8]-[C:4]-[C:5]-[I:6]>>[C:1]=[C:4].[I:3]-[C:2]([H:7])-[C:5]([H:8])-[I:6]')
        terminal='hydroxyl'
        sequence = 'AAAAAAAA'
        polymer = build.build_polymer(sequence = sequence,
                                        monomer_list = monomer_list,
                                        reaction = reaction,
                                        terminal = terminal)
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
        # Create a polymer system
        x = build.polymer_system(monomer_list=['O[C@H](C)C(=O)OI'], 
                    reaction = AllChem.ReactionFromSmarts('[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]'), 
                    length_target = 10, 
                    terminals = 'hydroxyl', 
                    num_chains = 5, 
                    copolymer=False)
        x.generate_conformers()
        # Calculate box components
        molecules, number_of_copies, topology, box_vectors, residual_monomer_actual = build.calculate_box_components(chains = x.chains,
                                                                                                                     monomers = x.monomers,
                                                                                                                    sequence=x.sequence)        

        # Check if the returned molecules is a list
        self.assertTrue(isinstance(molecules, list))
        self.assertEqual(len(molecules), 5)
        # Check if the returned number_of_copies is an integer
        self.assertTrue(isinstance(number_of_copies, list))
        # Check if the returned topology is a Topology object
        self.assertTrue(isinstance(topology, Topology))
        # Check if the returned box_vectors has the correct shape
        self.assertEqual(box_vectors.shape, (3, 3))   
        # Check if the returned residual_monomer_actual is a float
        self.assertTrue(isinstance(residual_monomer_actual, float))
        # Check if the returned residual_monomer_actual is as expected
        self.assertTrue(residual_monomer_actual==0.0)  
        x = build.polymer_system(monomer_list=['O[C@H](C)C(=O)OI'], 
                    reaction = AllChem.ReactionFromSmarts('[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]'), 
                    length_target = 10, 
                    terminals = 'hydroxyl', 
                    num_chains = 100, 
                    copolymer=False)
        x.generate_conformers()   

        # Calculate box components - test case with residual monomer
        molecules, number_of_copies, topology, box_vectors, residual_monomer_actual = build.calculate_box_components(chains = x.chains,
                                                                                                                     monomers = x.monomers,
                                                                                                                    sequence=x.sequence,
                                                                                                                    residual_monomer=0.5)        

        # Check if the returned molecules is a list
        self.assertTrue(isinstance(molecules, list))
        self.assertEqual(len(molecules), 5)
        # Check if the returned number_of_copies is an integer
        self.assertTrue(isinstance(number_of_copies, list))
        # Check if the returned topology is a Topology object
        self.assertTrue(isinstance(topology, Topology))
        # Check if the returned box_vectors has the correct shape
        self.assertEqual(box_vectors.shape, (3, 3))   
        # Check if the returned residual_monomer_actual is a float
        self.assertTrue(isinstance(residual_monomer_actual, float))
        # Check if the returned residual_monomer_actual is as expected
        self.assertAlmostEqual(round(residual_monomer_actual,3), 0.5, places=4)

class TestPolymerSystem(unittest.TestCase):
    def test_init(self):

        x = build.polymer_system(monomer_list=['O[C@H](C)C(=O)O[I]'], 
                    reaction = AllChem.ReactionFromSmarts('[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]'), 
                    length_target = 10, 
                    terminals = 'hydroxyl', 
                    num_chains = 5, 
                    copolymer=False)

        self.assertTrue(len(x.chains)==5)
        self.assertTrue(9 <= round(x.max_length)<= 11)
        self.assertTrue(1.5 <= round(x.PDI) <= 5)
        self.assertTrue(x.num_chains == 5)
        self.assertIsNotNone(x.monomers)
        x.generate_conformers()
        self.assertTrue(len(x.chains[0].conformers[0])==len(x.chains[0].atoms))
        #Test all charge methods
        x.charge_system('espaloma')
        self.assertTrue(len(x.chains[0].partial_charges)==len(x.chains[0].atoms))
        x.charge_system('NAGL')
        self.assertTrue(len(x.chains[0].partial_charges)==len(x.chains[0].atoms))

        #Check all chains have a name
        for chain in x.chains:
            assert chain.name is not None
        #Test case - Copolymer with 5% acceptance margin
        x = build.polymer_system(monomer_list=['O[C@H](C)C(=O)O[I]','OCC(=O)O[I]'], 
                                reaction = AllChem.ReactionFromSmarts('[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]'),
                                length_target=10,
                                num_chains = 5,
                                blockiness_target=1.0,
                                perc_A_target=50, 
                                copolymer=True,
                                acceptance=5)
        self.assertTrue(len(x.chains)==5)
        self.assertTrue(9 <= round(x.max_length)<= 11)
        self.assertTrue(47.5 <= x.A_actual <= 52.5)
        self.assertTrue(0.95 <= x.mean_blockiness <= 1.05)




# Run

if __name__ == '__main__':
    unittest.main()
