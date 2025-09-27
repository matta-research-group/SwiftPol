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

RDLogger.DisableLog("rdApp.*")


# Define all test cases
class TestBuildPolymer(unittest.TestCase):
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")

    def test_build_polymer(self):
        sequence = "AABBAABB"
        monomer_list = ["OC(=O)COI", "C[C@@H](C(=O)[OH])OI"]
        reaction = AllChem.ReactionFromSmarts(
            "[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]"
        )
        # Test the function
        polymer = build.build_polymer(
            sequence=sequence,
            monomer_list=monomer_list,
            reaction=reaction,
            terminal="hydrogen",
        )
        self.assertIsNotNone(polymer)
        # Test with an invalid sequence
        with self.assertRaises(IndexError):
            polymer = build.build_polymer("invalid", monomer_list, reaction)
        # Test with an invalid monomer list
        with self.assertRaises(IndexError):
            polymer = build.build_polymer(sequence, ["invalid"], reaction)
        # Test with an invalid reaction
        with self.assertRaises(AttributeError):
            polymer = build.build_polymer(sequence, monomer_list, "invalid")
        # Test with no terminal adjustment
        polymer = build.build_polymer(
            sequence="AAAAABBBBB",
            monomer_list=["O[C@H](C)C(=O)O[I]", "OCC(=O)O[I]"],
            reaction=reaction,
        )
        self.assertIsNotNone(polymer)
        for atom in polymer.GetAtoms():
            assert atom.GetPDBResidueInfo() is not None

        # Test ester terminals
        sequence = "AABBAABB"
        monomer_list = ["OC(=O)COI", "C[C@@H](C(=O)[OH])OI"]
        reaction = AllChem.ReactionFromSmarts(
            "[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]"
        )
        # Test the function
        polymer = build.build_polymer(
            sequence=sequence,
            monomer_list=monomer_list,
            reaction=reaction,
            terminal="ester",
        )
        self.assertIsNotNone(polymer)
        for atom in polymer.GetAtoms():
            assert atom.GetPDBResidueInfo() is not None

        # Test carboxyl terminals
        sequence = "AABBAABB"
        monomer_list = ["OC(=O)COI", "C[C@@H](C(=O)[OH])OI"]
        reaction = AllChem.ReactionFromSmarts(
            "[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]"
        )
        # Test the function
        polymer = build.build_polymer(
            sequence=sequence,
            monomer_list=monomer_list,
            reaction=reaction,
            terminal="carboxyl",
        )
        self.assertIsNotNone(polymer)
        for atom in polymer.GetAtoms():
            assert atom.GetPDBResidueInfo() is not None

        # Test SMILES terminals
        sequence = "AABBAABB"
        monomer_list = ["OC(=O)COI", "C[C@@H](C(=O)[OH])OI"]
        reaction = AllChem.ReactionFromSmarts(
            "[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]"
        )
        # Test the function
        polymer = build.build_polymer(
            sequence=sequence,
            monomer_list=monomer_list,
            reaction=reaction,
            terminal="CCCO[H]",
        )
        self.assertIsNotNone(polymer)

        # Test with cellulose
        monomer_list = [
            "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)OI)O)O)OI)O"
        ]  # Glucose with iodinated 1,4 positions
        reaction = AllChem.ReactionFromSmarts(
            "[I:1][O:2].[I:3][O:4][C:5]>>[C:5][O:2].[I:3][I:1].[O:4]"
        )
        sequence = "AAA"
        polymer = build.build_polymer(
            sequence=sequence,
            monomer_list=monomer_list,
            reaction=reaction,
            terminal="hydrogen",
        )
        self.assertIsNotNone(polymer)

        # Test with chitin
        monomer_list = ["C([C@@H]1[C@H]([C@@H]([C@H](C(O1)OI)N)O)OI)O"]
        reaction = AllChem.ReactionFromSmarts(
            "[I:1][O:2].[I:3][O:4][C:5]>>[C:5][O:2].[I:3][I:1].[O:4]"
        )
        terminal = "hydrogen"
        sequence = "AAA"
        polymer = build.build_polymer(
            sequence=sequence,
            monomer_list=monomer_list,
            reaction=reaction,
            terminal=terminal,
        )
        self.assertIsNotNone(polymer)

        # Test with polyethylene
        monomer_list = ["ICCI"]
        reaction = AllChem.ReactionFromSmarts(
            "[C:1][C:2][I:3].[C:4][C:5][I:6]>>[C:1][C:2][C:5][C:4].[I:3][I:6]"
        )
        terminal = "hydrogen"
        sequence = "AAAAAAAA"
        polymer = build.build_polymer(
            sequence=sequence,
            monomer_list=monomer_list,
            reaction=reaction,
            terminal=terminal,
        )
        self.assertIsNotNone(polymer)

        # Test with polythiophene
        monomer_list = ["S1C(CCI)=CC=C(CCI)1"]
        reaction = AllChem.ReactionFromSmarts(
            "[H:7]-[C:1]-[C:2]-[I:3].[H:8]-[C:4]-[C:5]-[I:6]>>[C:1]=[C:4].[I:3]-[C:2]([H:7])-[C:5]([H:8])-[I:6]"
        )
        terminal = "hydrogen"
        sequence = "AAAAAAAA"
        polymer = build.build_polymer(
            sequence=sequence,
            monomer_list=monomer_list,
            reaction=reaction,
            terminal=terminal,
        )
        self.assertIsNotNone(polymer)


class TestBuildLinearCopolymer(unittest.TestCase):
    def test_build_linear_copolymer(self):
        sequence = "ABAB"
        monomer_a_smiles = "OC(=O)CO"
        monomer_b_smiles = "C[C@@H](C(=O)[OH])O"
        reaction = AllChem.ReactionFromSmarts(
            "[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]"
        )
        polymer, A_ratio, B_ratio = build.build_linear_copolymer(
            sequence, monomer_a_smiles, monomer_b_smiles, reaction
        )

        # Check if the returned polymer is a RDKit Mol object
        self.assertTrue(isinstance(polymer, Chem.rdchem.Mol))

        # Check if the returned ratios are as expected
        self.assertEqual(A_ratio, 50.0)
        self.assertEqual(B_ratio, 50.0)


class TestPDI(unittest.TestCase):
    def test_PDI(self):
        # Create some RDKit molecule objects
        mol1 = Chem.MolFromSmiles("CC(=O)OC")
        mol2 = Chem.MolFromSmiles("CC(=O)OC")
        mol3 = Chem.MolFromSmiles("CC(=O)OC")

        chains = [mol1, mol2, mol3]

        pdi, mn, mw = build.PDI(chains)

        # Check if the returned PDI, Mn, and Mw are floats
        self.assertTrue(isinstance(pdi, float))
        self.assertTrue(isinstance(mn, float))
        self.assertTrue(isinstance(mw, float))

        # Check if the returned PDI, Mn, and Mw are as expected
        self.assertAlmostEqual(pdi, 1.00, places=2)
        self.assertAlmostEqual(mn, 74.0, places=1)
        self.assertAlmostEqual(mw, 74.0, places=1)


class TestBlockinessGen(unittest.TestCase):
    def test_blockiness_Gen(self):
        sequence = "ABABABAB"
        blockiness, block_length_A, block_length_B = build.blockiness_gen(sequence)

        # Check if the returned blockiness is float
        self.assertTrue(isinstance(blockiness, float))

        # Check if the returned blockiness, block_length_G, and block_length_L are as expected
        self.assertAlmostEqual(blockiness, 0.0, places=2)
        self.assertAlmostEqual(block_length_A, 1.0, places=2)
        self.assertAlmostEqual(block_length_B, 1.0, places=2)

        # Test case - PLA sequence
        sequence_A = "AAAAAAAA"
        blockiness_A, block_length_A_2, block_length_B_2 = build.blockiness_gen(
            sequence_A
        )
        self.assertTrue(isinstance(blockiness_A, str))

        # Test case - blockinesss wrt 'B', non-zero
        sequence = "AABBAABBAABB"
        blockiness, block_length_A, block_length_B = build.blockiness_gen(sequence, "B")
        self.assertAlmostEqual(blockiness, 1.0, places=2)

        # Test case - blockinesss wrt 'A', non-zero
        sequence = "AABBAABBAABB"
        blockiness, block_length_A, block_length_B = build.blockiness_gen(sequence, "A")
        self.assertAlmostEqual(blockiness, 1.5, places=2)


# Test calculate box components
class TestCalculateBoxComponents(unittest.TestCase):
    def test_calculate_box_components(self):
        from openff.units import unit

        # Create a polymer system
        x = build.polymer_system(
            monomer_list=["O[C@H](C)C(=O)O[I]", "OCC(=O)O[I]"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=30,
            terminals="hydrogen",
            num_chains=1,
            perc_A_target=75,
            blockiness_target=[1.0, "A"],
            copolymer=True,
            acceptance=1,
        )
        x.generate_conformers()

        molecules, number_of_copies, topology, box_vectors, residual_monomer_actual = (
            build.calculate_box_components(
                chains=x.chains,
                monomers=x.monomers,
                sequence=x.sequence,
                salt_concentration=0 * unit.mole / unit.liter,
                residual_monomer=1.5,
            )
        )

        # Check if the returned molecules is a list
        self.assertTrue(isinstance(molecules, list))
        self.assertEqual(len(molecules), 5)
        print(number_of_copies[0])
        self.assertTrue(number_of_copies[0] > 1)
        # Check if the returned number_of_copies is a list
        self.assertTrue(isinstance(number_of_copies, list))
        # Check if the returned topology is a Topology object
        self.assertTrue(isinstance(topology, Topology))
        # Check if the returned box_vectors has the correct shape
        self.assertEqual(box_vectors.shape, (3, 3))
        # Check if the returned residual_monomer_actual is a float
        self.assertTrue(isinstance(residual_monomer_actual, float))
        # Check if the returned residual_monomer_actual is as expected
        # print('resid_mon = ', residual_monomer_actual)
        # self.assertTrue(1.20 <= residual_monomer_actual <= 1.80)

        # Calculate box components - test case without residual monomer
        molecules, number_of_copies, topology, box_vectors, residual_monomer_actual = (
            build.calculate_box_components(
                chains=x.chains,
                monomers=x.monomers,
                sequence=x.sequence,
                salt_concentration=0 * unit.mole / unit.liter,
                residual_monomer=0,
            )
        )
        # Check if the returned molecules is a list of 5 objects, and water is not included
        self.assertTrue(isinstance(molecules, list))
        self.assertEqual(len(molecules), 5)
        self.assertEqual(len(number_of_copies), 5)
        self.assertTrue(number_of_copies[0] > 1)
        print(number_of_copies)
        # Check if the returned number_of_copies is an integer
        self.assertTrue(isinstance(number_of_copies, list))
        # Check if the returned topology is a Topology object
        self.assertTrue(isinstance(topology, Topology))
        # Check if the returned box_vectors has the correct shape
        self.assertEqual(box_vectors.shape, (3, 3))
        # Check if the returned residual_monomer_actual is a float
        self.assertTrue(isinstance(residual_monomer_actual, float))
        # Check if the returned residual_monomer_actual is as expected
        self.assertTrue(residual_monomer_actual == 0.0)


class TestPolymerSystem(unittest.TestCase):
    def test_init(self):

        x = build.polymer_system(
            monomer_list=["O[C@H](C)C(=O)O[I]"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=30,
            terminals="hydrogen",
            num_chains=5,
            copolymer=False,
        )
        self.assertTrue(x.num_chains == 5)
        # self.assertTrue(20 <= round(x.length_average)<= 40)
        self.assertIsNotNone(x.monomers)
        x.generate_conformers()
        self.assertTrue(len(x.chains[0].conformers[0]) == len(x.chains[0].atoms))
        # Test all charge methods
        x.charge_system("espaloma")
        self.assertTrue(len(x.chains[0].partial_charges) == len(x.chains[0].atoms))
        x.charge_system("NAGL")
        self.assertTrue(len(x.chains[0].partial_charges) == len(x.chains[0].atoms))

        # Check all chains have a name
        for chain in x.chains:
            assert chain.name is not None
        # Test case - Copolymer with 5% acceptance margin
        x = build.polymer_system(
            monomer_list=["O[C@H](C)C(=O)O[I]", "OCC(=O)O[I]"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=30,
            num_chains=5,
            blockiness_target=[1.0, "B"],
            perc_A_target=50,
            copolymer=True,
            acceptance=5,
        )
        self.assertTrue(len(x.chains) == 5)
        # self.assertTrue(20 <= round(x.length_average)<= 40)
        self.assertTrue(47.5 <= x.A_actual <= 52.5)
        self.assertTrue(0.95 <= x.mean_blockiness <= 1.05)
        self.assertTrue(1.0 < x.PDI < 2.5)

        # Test case with stereoisomers
        x = build.polymer_system(
            monomer_list=["O[C@H](C)C(=O)O[I]", "OCC(=O)O[I]"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=50,
            num_chains=5,
            blockiness_target=[1.0, "A"],
            perc_A_target=50,
            copolymer=True,
            acceptance=5,
            stereoisomerism_input=["A", 0.5, "O[C@@H](C)C(=O)O[I]"],
        )
        self.assertTrue(len(x.chains) == 5)
        # self.assertTrue(40 <= round(x.length_average)<= 60)
        self.assertTrue(47.5 <= x.A_actual <= 52.5)
        self.assertTrue(0.95 <= x.mean_blockiness <= 1.05)
        chain = x.chain_rdkit[0]
        # Check isomerism has occurred
        chiral_centers_after = Chem.FindMolChiralCenters(chain, includeUnassigned=True)
        stereocentres = [i[1] for i in chiral_centers_after]
        self.assertTrue(len(set(stereocentres)) == 2)

        # Test metadata export
        x.export_to_csv("polymer_system_data.csv")
        self.assertTrue(os.path.isfile("polymer_system_data.csv"))
        os.remove("polymer_system_data.csv")

        # Test packmol packing
        y = build.polymer_system(
            monomer_list=["O[C@H](C)C(=O)O[I]"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=10,
            num_chains=1,
            copolymer=False,
            acceptance=10,
        )
        y.generate_conformers()
        solvated_y = y.pack_solvated_system()
        self.assertIsNotNone(solvated_y)

        # Test polyply output
        y.charge_system("NAGL")
        y.generate_conformers()
        files = y.generate_polyply_files()

        for i in files:
            self.assertTrue(os.path.isfile(i))
            os.remove(i)
        os.remove("swiftpol_output_pointenergy.mdp")

        # Test residual calculation
        sys = build.polymer_system(
            monomer_list=["O[C@H](C)C(=O)O[I]", "OCC(=O)O[I]"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=50,
            terminals="hydrogen",
            num_chains=50,
            perc_A_target=75,
            copolymer=True,
            acceptance=10,
            stereoisomerism_input=["A", 0.5, "O[C@@H](C)C(=O)O[I]"],
        )
        sys.charge_system("NAGL")
        (
            molecules,
            number_of_copies,
            residual_monomer_actual,
            residual_oligomer_actual,
        ) = sys.calculate_residuals(residual_monomer=5, residual_oligomer=15)
        self.assertTrue(residual_monomer_actual >= 4 and residual_monomer_actual <= 6)
        self.assertTrue(
            residual_oligomer_actual >= 12 and residual_oligomer_actual <= 18
        )
        self.assertIsNotNone(molecules)
        self.assertIsNotNone(number_of_copies)
        self.assertTrue(len(molecules) == len(number_of_copies))


class TestPolymerSystemFromPDI(unittest.TestCase):
    def test_init(self):
        # Test PDI=1.7
        sys = build.polymer_system_from_PDI(
            monomer_list=["OC(=O)COI"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=50,
            num_chains=50,
            PDI_target=1.7,
            copolymer=False,
            acceptance=10,
        )
        self.assertTrue(len(sys.chains) == 50)
        self.assertTrue(40 <= sys.length_average <= 60)
        self.assertTrue(1.5 <= sys.PDI <= 2.0)
        # Test PDI=1.0 (monodisperse)
        sys = build.polymer_system_from_PDI(
            monomer_list=["OC(=O)COI"],
            reaction="[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]",
            length_target=50,
            num_chains=20,
            PDI_target=1.0,
            copolymer=False,
            acceptance=10,
        )
        self.assertTrue(0.9 <= round(sys.PDI,1) <= 1.1)
        # Test copolymer - removed whilst PDI calculation in dev
        x = build.polymer_system_from_PDI(monomer_list=['O[C@H](C)C(=O)O[I]','OCC(=O)O[I]'],
                                reaction = '[C:1][O:2][H:3].[I:4][O:5][C:6]>>[C:1][O:2][C:6].[H:3][O:5][I:4]',
                                length_target=50,
                                num_chains = 20,
                                perc_A_target=50,
                                PDI_target=1.0,
                                blockiness_target=[1.0, 'B'],
                                copolymer=True,
                                acceptance=20)
        self.assertIsNotNone(x)
        self.assertTrue(0.9 <= round(sys.PDI,1) <= 1.1)


# Run

if __name__ == "__main__":
    unittest.main()
