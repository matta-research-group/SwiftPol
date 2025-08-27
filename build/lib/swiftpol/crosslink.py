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
    for mol in halogens:
        hydrogen = Chem.MolFromSmiles('[H]')
        #Extract metadata from halogen to replace?
        mol_no_halogens = Chem.ReplaceSubstructs(mol, Chem.MolFromSmarts(mol), hydrogen)[0]
        Chem.SanitizeMol(mol_no_halogens)

    #Add explicit hydrogens
    mol_with_h = Chem.AddHs(mol_no_halogens)
    return mol_with_h

def build_branched_polymer(starting_polymer, 
                              num_iterations=1, 
                              probability_of_crosslinked_addition=0.5, 
                              probability_of_linear_addition=0.5,
                              reaction_templates=[('[Cl:1]-[*:2]-[*:3](-[I:4])-[*:7].[Cl:8]-[*:9]-[*:10](-[I:11])-[*:12].[F:15].[F:16]>>[*:12]-[*:10](-[F:15])-[*:9]-[*:2]-[*:3](-[F:16])-[*:7].[Cl:1].[I:4].[Cl:8].[I:11]', (Chem.MolFromSmiles('F'), Chem.MolFromSmiles('F'))),
                                                  ('[*:2]-[*:4](-[F:5])-[*:6]-[*:7]-[*:8](-[F:9])-[*:10].[Cl:13]-[*:14]-[*:15](-[I:16])-[*:17].[Br:20]>>[*:2]-[*:4](-[Br:20])-[*:6]-[*:7]-[*C:8](-[*:10])-[*:14]-[*:15](-[I:16])-[*:17].[Cl:1].[F:9].[F:5]', (Chem.MolFromSmiles('Br'),)),
                                                  ('[*:2]-[*:4](-[Br:5])-[*:6]-[*:7]-[*:8](-[*:9]-[*:12]-[*:13](-[I:14])-[*:15]).[Cl:18]-[*:19]-[*:20](-[I:21])-[*:22]>>[*:2]-[*:4]-[*:19]-[*:20](-[I:21])-[*:22]-[*:6]-[*:7]-[*:8](-[*:9]-[*:12]-[*:13](-[I:14])-[*:15]).[Br:5].[Cl:18]', ()),
                                                  ('[*:2]-[*:4](-[F:5])-[*:6]-[*:7]-[*:8](-[F:9])-[*:10].[*:14]-[*:16](-[F:17])-[*:18]-[*:19]-[*:20](-[F:21])-[*:22]>>[*:2]-[*:4](-[*:16](-[*:14])(-[*:18]-[*:19]-[*:20](-[F:21])-[*:22]))-[*:6]-[*:7]-[*:8](-[F:9])-[*:10].[F:5].[F:17]', ()),]
                                                  ):
    """
    Builds a crosslinked polymer network by iteratively applying chemical reactions.


    This function builds a branched polymer by either adding linear chains 
    or crosslinking existing polymer chains. The user can specify the number of iterations 
    and the probabilities of crosslinking versus linear chain addition.

    Parameters
    ----------
    starting_polymer : rdkit.Chem.Mol
        The initial polymer molecule/monomer to start the crosslinking or chain addition process.
    num_iterations : int, optional
        The number of reaction iterations to perform (default is 1).
    probability_of_crosslinked_addition : float, optional
        The probability of performing a crosslinking reaction in each iteration (default is 0.5).
    probability_of_linear_addition : float, optional
        The probability of performing a linear chain addition reaction in each iteration (default is 0.5).

    Returns
    -------
    rdkit.Chem.Mol
        The final polymer molecule after the specified number of iterations, with halogens replaced by hydrogens.

    Raises
    ------
    AssertionError
        If the sum of `probability_of_crosslinked_addition` and `probability_of_linear_addition` is not equal to 1.

    Notes
    -----
    - The function uses predefined reaction templates (SMARTS) for crosslinking and chain addition.
    - Crosslinking requires at least two polymer chains to be present in the history.
    - The function randomly selects reactions and reactants based on the specified probabilities.
    - Halogens in the final polymer are replaced with hydrogens for chemical validity.

    Examples
    --------
    >>> from rdkit import Chem
    >>> starting_polymer = build.build_starting_polymer()
    >>> final_polymer = build_crosslinked_polymer(starting_polymer, num_iterations=5, 
    ...                                           probability_of_crosslinked_addition=0.6, 
    ...                                           probability_of_linear_addition=0.4)
    """
#### TO DO:
#### - rebuild this function with a simpler reaction template system.
#### - should be 2 types of templates - a dictionary for branched and linear?
#### - a target Mw for the polymer rather than iterations? - more intuitive. alternatively number of chains
#### - test reaction smarts first, and throw error if they don't work.
### - get rid of 'history' stuff
#### - increase residue number for each added chain/monomer

    from swiftpol.crosslink import replace_halogens_with_hydrogens
    # Check if the probabilities sum to 1
    if probability_of_crosslinked_addition + probability_of_linear_addition != 1:
        raise AssertionError("Probabilities must sum to 1.")

    #Convert SMARTS strings into RDKit reactions
    reactions = [(AllChem.ReactionFromSmarts(smarts), reactants) for smarts, reactants in reaction_templates]
        
    polymer_network = starting_polymer  #This is the main polymer we modify
    polymer_history = [polymer_network]  #Stores each iteration of the polymer network

    successful_reactions = 0
    while successful_reactions < num_iterations:
        #Choose whether to add a chain or crosslink
        action = random.choices(["add_chain", "crosslink"], weights=[probability_of_linear_addition, probability_of_crosslinked_addition], k=1)[0]

        if action == "crosslink":
            #Ensure crosslinking only happens when there's something to crosslink with
            if len(polymer_history) < 2:
                continue  #Skip this step and retry

            #Choose a previous polymer to crosslink with
            chosen_polymer = random.choice(polymer_history[:-1])  #Select from earlier polymer versions
            chosen_reaction, chosen_reactants = random.choices(reactions, weights=[0.0,0.3,0.35,0.35], k=1)[0]
            reactant_sequence = (polymer_network, chosen_polymer) + chosen_reactants  #Use current and chosen polymer

        else:  #Adding a chain
            chosen_reaction, chosen_reactants = reactions[0]
            reactant_sequence = (polymer_network, starting_polymer) + chosen_reactants  #Add chain to current polymer

        #Apply the reaction
        reactant_sets = chosen_reaction.RunReactants(reactant_sequence)

        if not reactant_sets:
            continue  #Retry without counting this attempt

        #Choose one of the successful products
        polymer_network = random.choice(reactant_sets)[0]
        polymer_history.append(polymer_network)  #Save the new polymer version
        successful_reactions += 1  #Count successful reactions

    return replace_halogens_with_hydrogens(polymer_history[-1])  #Return last polymer in the history


def assign_random_stereochemistry(mol):
    """
    Takes an RDKit Mol with undefined stereocenters and returns a copy
    with random stereochemistry assigned to all stereocenters.
    """
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
    from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions
    opts = StereoEnumerationOptions(
        onlyUnassigned = True,              #"Only look at stereocenters that don't have a defined configuration (i.e. ambiguous centers)."
        unique = True                       #"Don't generate redundant stereoisomers - give me distinct ones."
    )

    isomers = list(EnumerateStereoisomers(mol, options = opts))   
    #Runs a stereo enumeration algorithm ---> If there are 'n' undefined stereocenters, this could generate up to 2^n stereoisomers.
    #Prepares the list of all the options.
    
    #If any isomers were generated -> Randomly pick one. Else return the original molecule unchanged.
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
    rw = Chem.RWMol(mol)                       #Let's us add/remove atoms and bonds
    cl_sites = find_terminal_cl_sites(rw)

    pairs = list(zip(cl_sites[::2], cl_sites[1::2]))               #As above
    if max_pairs is not None:                                      #If max pairs is provided we override pairing
        pairs = filter_close_pairs(cl_sites, max_index_diff = 150)

    removed_indices = set()

    #Step 1: Add all bonds first
    for (cl1, c1), (cl2, c2) in pairs:                                 #Creates a bond between carbon atoms
        #Skip any atoms that are already deleted in this loop
        if cl1 in removed_indices or cl2 in removed_indices:
            continue
        try:
            rw.AddBond(c1, c2, Chem.BondType.SINGLE)
        except:
            print(f"Failed to add bond between atoms {c1} and {c2}")
        
    #Step 2: Remove Cl atoms after all bonds
    atoms_to_remove = sorted([cl for (cl1, _), (cl2, _) in pairs for cl in (cl1, cl2)], reverse = True)
    for cl in atoms_to_remove:
        if cl not in removed_indices:
            #Remove Cl atoms first
            rw.RemoveAtom(cl)
            removed_indices.add(cl)

    Chem.SanitizeMol(rw)
    return rw

def iterative_crosslink(mol, max_cycles = 10):
    cycles = 0
    history = [mol]
    while cycles < max_cycles:
        cl_sites = find_terminal_cl_sites(mol)
        if len(cl_sites) < 2:
            break
        mol = crosslink_cl_sites(mol)
        history.append(mol)
        cycles += 1
    return mol, history


def crosslink_polymer(branched_polymer):
    '''
    docstring
    '''
### remember to add imports
    return crosslinked_polymer

# crosslinked_polymer class?