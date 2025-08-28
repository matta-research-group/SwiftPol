import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
from swiftpol import build
from swiftpol.crosslink import replace_halogens_with_hydrogens, validate_linear_reaction
from swiftpol import crosslink


class TestCrosslink(unittest.TestCase):
    def test_replace_halogens_with_hydrogens(self):
        """Test that halogens are replaced with hydrogens in the molecule."""
        # Input molecule with halogens
        mol = Chem.MolFromSmiles("ClCCBr")
        mol = Chem.AddHs(mol)
        # Add residue info to atoms for residue copying testing
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName('TEST')
        info.SetResidueNumber(1)
        info.SetChainId('A')
        [atom.SetMonomerInfo(info) for atom in mol.GetAtoms()]
        self.assertIsNotNone(mol, "Failed to create molecule from SMILES.")

        # Replace halogens with hydrogens
        modified_mol = replace_halogens_with_hydrogens(mol)

        # Check that no halogens remain
        halogens = ["F", "Cl", "Br", "I"]
        for halogen in halogens:
            halogen_smarts = Chem.MolFromSmarts(halogen)
            self.assertFalse(
                modified_mol.HasSubstructMatch(halogen_smarts),
                f"Halogen {halogen} was not removed.",
            )

        # Check that hydrogens were added
        num_hydrogens = sum(1 for atom in modified_mol.GetAtoms() if atom.GetSymbol() == "H")
        self.assertGreater(num_hydrogens, 0, "No hydrogens were added to the molecule.")

    def test_replace_halogens_with_hydrogens_no_halogens(self):
        """Test that the function works correctly when no halogens are present."""
        # Input molecule without halogens
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName('TEST')
        info.SetResidueNumber(1)
        info.SetChainId('A')
        [atom.SetMonomerInfo(info) for atom in mol.GetAtoms()]
        self.assertIsNotNone(mol, "Failed to create molecule from SMILES.")

        # Replace halogens with hydrogens
        modified_mol = replace_halogens_with_hydrogens(mol)

        # Check that the molecule remains unchanged
        self.assertEqual(
            Chem.MolToSmiles(modified_mol),
            Chem.MolToSmiles(mol),
            "Molecule was modified even though no halogens were present.",
        )

    def test_validate_linear_reaction_valid(self):
        """Test that validate_linear_reaction works with valid reactions."""
        # Create a simple starting polymer
        monomer = ["IC=CC(=O)OI", "IOCCOI"] # PEGDA monomers with iodinated polymerization points
        reaction = AllChem.ReactionFromSmarts('[C:1]-[O:2]-[I:3].[C:4]-[O:5]-[I:6]>>[C:1]-[O:2]-[C:4].[I:3]-[I:6].[O:5]')
        sequence = 'ABBBBBBBBBA' # SwiftPol can easily build irregular sequence motifs
        starting_polymer = build.build_polymer(sequence = sequence,
                                        monomer_list = monomer,
                                        reaction = reaction,
                                        terminal = 'hydrogen',
                                        chainID='A')

        self.assertIsNotNone(starting_polymer, "Failed to create starting polymer.")

        # Create a linear activation reaction
        linear_activate = AllChem.ReactionFromSmarts('[C:1]=[C:2].[Cl:3].[I:4]>>[C:1](-[I:4])-[C:2]-[Cl:3]')
        self.assertIsNotNone(linear_activate, "Failed to create linear activation reaction.")

        # Create reactants for the activation reaction
        linear_activate_reactants = [Chem.MolFromSmiles("Cl"), Chem.MolFromSmiles("I")]
        self.assertTrue(all(linear_activate_reactants), "Failed to create activation reactants.")

        # Create a linear polymerization reaction
        linear_react = AllChem.ReactionFromSmarts('[C:1][I:2].[C:3][Cl:4]>>[C:1][C:3].[I:2].[Cl:4]')
        self.assertIsNotNone(linear_react, "Failed to create linear polymerization reaction.")

        # Validate the reactions
        try:
            validate_linear_reaction(starting_polymer, linear_activate, linear_activate_reactants, linear_react)
        except ValueError as e:
            self.fail(f"validate_linear_reaction raised ValueError unexpectedly: {e}")

    def test_validate_linear_reaction_invalid(self):
        """Test that validate_linear_reaction raises an error for invalid reactions."""
        # Create a simple starting polymer
        starting_polymer = Chem.MolFromSmiles("CCO")
        self.assertIsNotNone(starting_polymer, "Failed to create starting polymer.")

        # Create an invalid linear activation reaction
        linear_activate = AllChem.ReactionFromSmarts("[C:1][O:2]>>[C:1][N:2]")  # Invalid reaction
        self.assertIsNotNone(linear_activate, "Failed to create linear activation reaction.")

        # Create reactants for the activation reaction
        linear_activate_reactants = [Chem.MolFromSmiles("CCO")]
        self.assertTrue(all(linear_activate_reactants), "Failed to create activation reactants.")

        # Create a linear polymerization reaction
        linear_react = AllChem.ReactionFromSmarts("[C:1][O:2].[C:3][O:4]>>[C:1][O:2][C:3][O:4]")
        self.assertIsNotNone(linear_react, "Failed to create linear polymerization reaction.")

        # Validate the reactions and expect a ValueError
        with self.assertRaises(ValueError, msg="validate_linear_reaction did not raise ValueError for invalid reactions."):
            validate_linear_reaction(starting_polymer, linear_activate, linear_activate_reactants, linear_react)


    def test_build_branched_polymer(self):
        monomer = ["IC=CC(=O)OI", "IOCCOI"] # PEGDA monomers with iodinated polymerization points
        reaction = AllChem.ReactionFromSmarts('[C:1]-[O:2]-[I:3].[C:4]-[O:5]-[I:6]>>[C:1]-[O:2]-[C:4].[I:3]-[I:6].[O:5]')
        sequence = 'ABBBBBBBBBA' # SwiftPol can easily build irregular sequence motifs
        polymer = build.build_polymer(sequence = sequence,
                                        monomer_list = monomer,
                                        reaction = reaction,
                                        terminal = 'hydrogen',
                                        chainID='A')
        
        reaction_templates = {'linear_activate' : ['[C:1]=[C:2].[Cl:3].[I:4]>>[C:1](-[I:4])-[C:2]-[Cl:3]', 'Cl', 'I'],
                      'linear_react' : ['[C:1][I:2].[C:3][Cl:4]>>[C:1][C:3].[I:2].[Cl:4]'],
                      'branched_activate' : ['[O:1][C:2][C:3](-[C:4][Cl:5])[C:6][C:7][C:8].[Br:9]>>[O:1][C:2][C:3](-[C:4][Br:9])[C:6][C:7][C:8].[Cl:5]', 'Br'],
                      'branched_react' : ['[C:1][Br:2].[C:3][Cl:4]>>[C:1][C:3].[Br:2].[Cl:4]']}
        
        branched_polymer = crosslink.build_branched_polymer(starting_polymer = polymer,
                                                            reaction_templates = reaction_templates,
                                                            num_iterations=4, 
                                                            probability_of_branched_addition=0.5, 
                                                            probability_of_linear_addition=0.5)
        
        self.assertIsNotNone(branched_polymer, "Failed to create branched polymer.")

        branched_polymer_mw = crosslink.build_branched_polymer(starting_polymer = polymer,
                                                                reaction_templates = reaction_templates,
                                                                target_mw=10000, 
                                                                probability_of_branched_addition=0.5, 
                                                                probability_of_linear_addition=0.5)

        self.assertIsNotNone(branched_polymer_mw, "Failed to create branched polymer with target MW.")

        def test_build_crosslinked_polymer(self):
            polymer = build.build_polymer(sequence = sequence,
                                        monomer_list = monomer,
                                        reaction = reaction,
                                        terminal = 'hydrogen',
                                        chainID='A')
        
        
            reaction_templates = {'linear_activate' : ['[C:1]=[C:2].[Cl:3].[I:4]>>[C:1](-[I:4])-[C:2]-[Cl:3]', 'Cl', 'I'],
                                'linear_react' : ['[C:1][I:2].[C:3][Cl:4]>>[C:1][C:3].[I:2].[Cl:4]'],
                                'branched_activate' : ['[O:1][C:2][C:3](-[C:4][Cl:5])[C:6][C:7][C:8].[Br:9]>>[O:1][C:2][C:3](-[C:4][Br:9])[C:6][C:7][C:8].[Cl:5]', 'Br'],
                                'branched_react' : ['[C:1][Br:2].[C:3][Cl:4]>>[C:1][C:3].[Br:2].[Cl:4]']}
        
            branched_polymer = crosslink.build_branched_polymer(starting_polymer = polymer,
                                                                reaction_templates = reaction_templates,
                                                                num_iterations=4, 
                                                                probability_of_branched_addition=0.5, 
                                                                probability_of_linear_addition=0.5)

            crosslinked_network = crosslink.crosslink_polymer(branched_mw, percentage_to_crosslink=80)
            self.assertIsNotNone(crosslinked_network, "Failed to create crosslinked polymer network.")



if __name__ == "__main__":
    unittest.main()