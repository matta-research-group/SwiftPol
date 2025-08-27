##Functions for building crosslinked polymers

#Init
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import random
import numpy as np
try:
    from openeye import oechem
    oechem_imported = True
except:
    import warnings
    warnings.warn("OpenEye is not installed. You will not be able to use OpenEye Toolkits for conformer generation.")
    oechem_imported = False



def replace_halogens_with_hydrogens(mol):
    """
    Replaces all halogen atoms (F, Cl, Br, I) in an RDKit molecule with hydrogen atoms.

    This function iterates through the atoms in the given molecule, identifies halogen atoms 
    (fluorine, chlorine, bromine, iodine etc.), removes them, and replaces them with explicit hydrogens. 

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        The input RDKit molecule from which halogens will be replaced with hydrogens.

    Returns
    -------
    rdkit.Chem.Mol
        A new RDKit molecule with all halogens replaced by explicit hydrogens.

    """
    halogens = ['F', 'Cl', 'Br', 'I', 'At', 'Ts']
    for halogen in halogens:
        try:
            matches = mol.GetSubstructMatches(Chem.MolFromSmarts(halogen))
        except:
            matches = []
        if len(matches) > 0:
            for match in matches:
                #extract pdb info of neighbouring atom and assign to hydrogen
                atom = mol.GetAtomWithIdx(match[0])
                bonded_atom = atom.GetNeighbors()
                info = bonded_atom[0].GetMonomerInfo()
                info_new = Chem.AtomPDBResidueInfo()
                info_new.SetResidueName(info.GetResidueName())
                info_new.SetResidueNumber(info.GetResidueNumber())
                info_new.SetChainId(info.GetChainId())
                hydrogen = Chem.MolFromSmiles('[H]')
                [atom.SetMonomerInfo(info_new)  for  atom  in  hydrogen.GetAtoms()]
                #Replace halogen with hydrogen
                mw = Chem.RWMol(mol)
                mw.ReplaceAtom(match[0],list(hydrogen.GetAtoms())[0])
                Chem.SanitizeMol(mw)
                mw.CommitBatchEdit()
                mol = Chem.AddHs(mw)
        
    mol_no_halogens = mol
    Chem.SanitizeMol(mol_no_halogens)
    #Add explicit hydrogens
    mol_with_h = Chem.AddHs(mol_no_halogens)
    return mol_with_h

from rdkit import Chem

def validate_linear_reaction(starting_polymer, linear_activate, linear_activate_reactants, linear_react):
    try:
        # Test the linear activation reaction
        activated_pol_test = linear_activate.RunReactants([starting_polymer] + linear_activate_reactants)
        if not activated_pol_test:
            raise ValueError("Linear activation reaction produced no results.")
        activated_pol_test = activated_pol_test[0][0]
        Chem.SanitizeMol(activated_pol_test)
    except Exception as e:
        raise ValueError(
            "Linear activation reaction SMARTS is invalid or incompatible with the starting polymer or linear activation reactants. "
            "For support with constructing reaction SMARTS, raise an issue at https://github.com/matta-research-group/SwiftPol/issues"
        ) from e

    try:
        # Test the linear polymerization reaction
        linear_polymer_test = linear_react.RunReactants([activated_pol_test, activated_pol_test])
        if not linear_polymer_test:
            raise ValueError("Linear polymerization reaction produced no results.")
        linear_polymer_test = linear_polymer_test[0][0]
        Chem.SanitizeMol(linear_polymer_test)
    except Exception as e:
        raise ValueError(
            "Linear reaction SMARTS is invalid or incompatible with the starting polymer. "
            "For support with constructing reaction SMARTS, raise an issue at https://github.com/matta-research-group/SwiftPol/issues"
        ) from e
    
def validate_branched_reaction(starting_polymer, branched_activate, branched_activate_reactants, branched_react, linear_activate, linear_activate_reactants):
    try:
        # Test the branched activation reaction
        activated_pol_test = branched_activate.RunReactants([starting_polymer] + branched_activate_reactants)
        if not activated_pol_test:
            raise ValueError("Branched activation reaction produced no results.")
        activated_pol_test = activated_pol_test[0][0]
        Chem.SanitizeMol(activated_pol_test)
    except Exception as e:
        raise ValueError(
            "Branched activation reaction SMARTS is invalid or incompatible with the polymer network or branched activation reactants. "
            "For support with constructing reaction SMARTS, raise an issue at https://github.com/matta-research-group/SwiftPol/issues"
        ) from e

    try:
        # Test the branched polymerization reaction
        linear_activated = linear_activate.RunReactants([starting_polymer] + linear_activate_reactants)[0][0]
        Chem.SanitizeMol(linear_activated)
        branched_polymer_test = branched_react.RunReactants([activated_pol_test, linear_activated])
        if not branched_polymer_test:
            raise ValueError("Branched polymerization reaction produced no results.")
        branched_polymer_test = branched_polymer_test[0][0]
        Chem.SanitizeMol(branched_polymer_test)
    except Exception as e:
        raise ValueError(
            "Branched reaction SMARTS is invalid or incompatible with the polymer network. "
            "For support with constructing reaction SMARTS, raise an issue at https://github.com/matta-research-group/SwiftPol/issues"
        ) from e

def iterative_chainID_update(polymer, last_polymer_added):
    import string
    for i in last_polymer_added.GetAtoms():
        index = i.GetIdx()
        info = i.GetMonomerInfo()
        if info is not None:
            info_new = Chem.AtomPDBResidueInfo()
            info_new.SetResidueName(f'{int(info.GetResidueName()[0])+1}{info.GetResidueName()[1:]}')
            info_new.SetResidueNumber(info.GetResidueNumber())
            info_new.SetChainId(string.ascii_uppercase[(string.ascii_uppercase.index(info.GetChainId()) + 1) % 26])
            polymer.GetAtomWithIdx(index).SetMonomerInfo(info_new)
        else:
            print(f"Warning: Atom {index} ({i.GetSymbol()}) does not have MonomerInfo. Skipping update for this atom.")
    return polymer

def build_branched_polymer(starting_polymer, 
                           reaction_templates,
                           num_iterations=None,
                           target_mw=None, 
                           probability_of_branched_addition=0.5, 
                           probability_of_linear_addition=0.5):
    def build_branched_polymer(starting_polymer, 
                               reaction_templates,
                               num_iterations=None,
                               target_mw=None, 
                               probability_of_branched_addition=0.5, 
                               probability_of_linear_addition=0.5):
        """
        Builds a branched polymer.
    
        Parameters:
        -----------
        starting_polymer : rdkit.Chem.Mol
            The initial polymer molecule to which reactions will be applied.
    
        reaction_templates : dict
            A dictionary containing reaction templates for both linear and branched additions.
            Keys should include:
                - "branched_activate": Activation reaction for branching.
                - "branched_react": Reaction for branching.
                - "linear_activate": Activation reaction for linear addition.
                - "linear_react": Reaction for linear addition.
    
        num_iterations : int, optional
            The number of iterations to perform. If specified, the process will stop after this many
            successful reactions, regardless of the target molecular weight.
    
        target_mw : float, optional
            The target molecular weight for the polymer. If specified, the process will stop once
            the molecular weight of the polymer reaches or exceeds this value.
    
        probability_of_branched_addition : float, optional, default=0.5
            The probability of performing a branched addition at each iteration.
    
        probability_of_linear_addition : float, optional, default=0.5
            The probability of performing a linear addition at each iteration.
    
        Returns:
        --------
        polymer_network : rdkit.Chem.Mol
            The resulting polymer molecule after the specified number of iterations or upon reaching
            the target molecular weight.
            Random stereochemsitry assigned to each stereocenter.
    
        Notes:
        ------
        - The function alternates between linear and branched additions based on the specified probabilities.
        - If both `num_iterations` and `target_mw` are specified, the process will stop as soon as either
          condition is met.
        - The reaction templates must be compatible with the starting polymer and follow RDKit's reaction
          SMARTS format.
    
        Raises:
        -------
        ValueError
            If the reaction templates are invalid or incompatible with the starting polymer.
        AssertionError
            If the sum of probabilities does not equal 1.
        ValueError
            If neither `num_iterations` nor `target_mw` is specified.
        """
    from swiftpol.crosslink import validate_linear_reaction, validate_branched_reaction, iterative_chainID_update, assign_random_stereochemistry
    # Check if the probabilities sum to 1
    if probability_of_branched_addition + probability_of_linear_addition != 1:
        raise AssertionError("Probabilities must sum to 1.")
    # Check that there is either a target Mw or a number of iterations
    if target_mw is None and num_iterations is None:
        raise ValueError("Either target_mw or num_iterations must be specified.")
    if target_mw is not None and num_iterations is not None:
        raise ValueError("Specify either target_mw or num_iterations, not both.")

    #Convert SMARTS strings into RDKit reactions
    rxns = list(reaction_templates.keys()) 
    linear_activate = AllChem.ReactionFromSmarts(reaction_templates[rxns[0]][0])
    linear_activate_reactants = [Chem.MolFromSmiles(smiles) for smiles in reaction_templates[rxns[0]][1:]]
    linear_react = AllChem.ReactionFromSmarts(reaction_templates[rxns[1]][0])
    validate_linear_reaction(starting_polymer, linear_activate, linear_activate_reactants, linear_react) #validate linear reactions

    branched_activate = AllChem.ReactionFromSmarts(reaction_templates[rxns[2]][0])
    branched_activate_reactants = [Chem.MolFromSmiles(smiles) for smiles in reaction_templates[rxns[2]][1:]]
    branched_react = AllChem.ReactionFromSmarts(reaction_templates[rxns[3]][0])
    validate_branched_reaction(starting_polymer, branched_activate, branched_activate_reactants, branched_react, linear_activate, linear_activate_reactants) #validate branched reactions
    
    starting_polymer_activated = linear_activate.RunReactants([starting_polymer] + linear_activate_reactants)[0][0]
    Chem.SanitizeMol(starting_polymer_activated)
    starting_polymer_activated = linear_activate.RunReactants([starting_polymer_activated] + linear_activate_reactants)[0][0]
    Chem.SanitizeMol(starting_polymer_activated)
    Chem.AddHs(starting_polymer_activated)
    polymer_network = Chem.Mol(starting_polymer_activated)
    last_polymer_added = Chem.Mol(starting_polymer)
    if num_iterations is not None:
        successful_reactions = 0
        while successful_reactions < num_iterations:
            # Choose whether to complete a linear addition or crosslink
            action = random.choices(["linear", "crosslink"], weights=[probability_of_linear_addition, probability_of_branched_addition], k=1)[0]
            if action == "crosslink":
                activated_for_branching_results = branched_activate.RunReactants([polymer_network] + branched_activate_reactants)
                if activated_for_branching_results:  
                    n = random.randint(0, len(activated_for_branching_results) - 1)
                    activated_for_branching = activated_for_branching_results[n][0]
                    Chem.SanitizeMol(activated_for_branching)
                    starting_polymer_updated = iterative_chainID_update(starting_polymer_activated, last_polymer_added) #Update ChainID
                    branched_polymer_results = branched_react.RunReactants([activated_for_branching, starting_polymer_updated])
                    if branched_polymer_results:
                        last_polymer_added = Chem.Mol(starting_polymer_updated)
                        branched_polymer = branched_polymer_results[0][0]
                        Chem.SanitizeMol(branched_polymer)
                        polymer_network = branched_polymer
                        successful_reactions += 1

            elif action == "linear":
                starting_polymer_updated = iterative_chainID_update(starting_polymer_activated, last_polymer_added) #Update ChainID
                linear_polymer_results = linear_react.RunReactants([polymer_network, starting_polymer_updated])
                if linear_polymer_results: 
                    last_polymer_added = Chem.Mol(starting_polymer_updated)
                    n = random.randint(0, len(linear_polymer_results) - 1)
                    linear_polymer = linear_polymer_results[n][0]
                    Chem.SanitizeMol(linear_polymer)
                    polymer_network = linear_polymer
                    successful_reactions += 1

    elif target_mw is not None:
        mol_weight_rolling = Chem.Descriptors.MolWt(polymer_network)
        while mol_weight_rolling < target_mw:
            # Choose whether to add a chain or crosslink
            action = random.choices(["linear", "crosslink"], weights=[probability_of_linear_addition, probability_of_branched_addition], k=1)[0]
    
            if action == "crosslink":
                activated_for_branching_results = branched_activate.RunReactants([polymer_network] + branched_activate_reactants)
                if activated_for_branching_results:  
                    n = random.randint(0, len(activated_for_branching_results) - 1)
                    activated_for_branching = activated_for_branching_results[n][0]
                    Chem.SanitizeMol(activated_for_branching)
                    branched_polymer_results = branched_react.RunReactants([activated_for_branching, polymer_network])
                    if branched_polymer_results: 
                        branched_polymer = branched_polymer_results[0][0]
                        Chem.SanitizeMol(branched_polymer)
                        polymer_network = branched_polymer
    
            elif action == "linear":
                linear_polymer_results = linear_react.RunReactants([polymer_network, starting_polymer])
                if linear_polymer_results: 
                    n = random.randint(0, len(linear_polymer_results) - 1)
                    linear_polymer = linear_polymer_results[n][0]
                    Chem.SanitizeMol(linear_polymer)
                    polymer_network = linear_polymer
            mol_weight_rolling = Chem.Descriptors.MolWt(polymer_network)

    return assign_random_stereochemistry(polymer_network) 


def assign_random_stereochemistry(mol):
    """
    Takes an RDKit Mol with undefined stereocenters and returns a copy
    with random stereochemistry assigned to all stereocenters.
    """
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
    from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions
    opts = StereoEnumerationOptions(
        onlyUnassigned = True,              
        unique = True                      
    )

    isomers = list(EnumerateStereoisomers(mol, options = opts))   
    
    return random.choice(isomers) if isomers else mol

###To do here:
# ---- logic works but needs to be streamlined
# --- One function to find leftover terminals/add halogens/
# ---- How does this relate to degree of crosslinking?
# Then clear halogens
# another function to create a network?

def find_terminal_cl_sites(mol):
    cl_sites = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "Cl":
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == "C":
                cl_sites.append((atom.GetIdx(), neighbor.GetIdx()))
    return cl_sites


def filter_close_pairs(cl_sites, max_index_diff = 150):
    return [((cl1, c1), (cl2, c2))
            for (cl1, c1), (cl2, c2) in zip(cl_sites[::2], cl_sites[1::2])     #Pairs up adjacent elements [(0, 1), (2, 3)] -> [((0,1), (2,3))]
            if abs(c1 - c2) <= max_index_diff]                                 #For each pair it checks the distance between the two atoms

def crosslink_cl_sites(mol, max_pairs = None):
    rw = Chem.RWMol(mol)        
    cl_sites = find_terminal_cl_sites(rw)

    pairs = list(zip(cl_sites[::2], cl_sites[1::2]))               #As above
    if max_pairs is not None:                                      #If max pairs is provided we override pairing
        pairs = filter_close_pairs(cl_sites, max_index_diff = 150)

    removed_indices = set()

    #Add all bonds first
    for (cl1, c1), (cl2, c2) in pairs:                                 #Creates a bond between carbon atoms
        #Skip any atoms that are already deleted in this loop
        if cl1 in removed_indices or cl2 in removed_indices:
            continue
        try:
            rw.AddBond(c1, c2, Chem.BondType.SINGLE)
        except:
            print(f"Failed to add bond between atoms {c1} and {c2}")
        
    #Remove Cl atoms after all bonds
    atoms_to_remove = sorted([cl for (cl1, _), (cl2, _) in pairs for cl in (cl1, cl2)], reverse = True)
    for cl in atoms_to_remove:
        if cl not in removed_indices:
            #Remove Cl atoms first
            rw.RemoveAtom(cl)
            removed_indices.add(cl)

    Chem.SanitizeMol(rw)
    rw.CommitBatchEdit()
    mol = Chem.AddHs(rw)
    return mol

def iterative_crosslink(mol, max_cycles = 10):
    from swiftpol.crosslink import find_terminal_cl_sites, crosslink_cl_sites
    cycles = 0
    while cycles < max_cycles:
        cl_sites = find_terminal_cl_sites(mol)
        if len(cl_sites) < 2:
            break
        mol = crosslink_cl_sites(mol)
        cycles += 1
    return mol


def crosslink_polymer(branched_polymer):
    '''
    docstring
    '''
### remember to add imports
    return crosslinked_polymer

# crosslinked_polymer class?