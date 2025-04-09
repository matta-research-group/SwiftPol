#Init
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import random
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
try:
    from openeye import oechem
    oechem_imported = True
except:
    import warnings
    warnings.warn("OpenEye is not installed. You will not be able to use OpenEye Toolkits for conformer generation.")

import openmm
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import get_data_file_path
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, AmberToolsToolkitWrapper
from openff.units import unit
from pandas import read_csv


from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box, _max_dist_between_points


from functools import reduce
from statistics import mean
from rdkit.Chem.Descriptors import ExactMolWt
from openff.interchange import Interchange


#Build polymer - generic
def build_polymer(sequence, monomer_list, reaction, terminal='hydroxyl', chain_num=1):
    """
    Constructs a polymer from a given sequence of monomers.

    Parameters
    ----------
    sequence : str
        A string representing the sequence of monomers (e.g., 'ABAB').
    monomer_list : list
        A list of SMILES strings representing the monomers.
    reaction : rdkit.Chem.rdChemReactions.ChemicalReaction
        An RDKit reaction object used to link monomers.
    terminal : str, optional
        The terminal group to be added to the polymer. Options are 'hydroxyl', 'carboxyl', or 'ester'.
        Default is 'hydroxyl'.
    chain_number : int, optional
        The number of polymer chains to construct. Default is 1. Input used for ensemble build.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        The constructed polymer as an RDKit molecule object.
    """
    from rdkit import RDLogger 
    RDLogger.DisableLog('rdApp.*')   
    monomers = {}
    for x in sorted(list(set(sequence))):
        ind = sorted(list(set(sequence))).index(x)
        monomers[x] = monomer_list[ind]
    hits = Chem.MolFromSmiles(monomers[sequence[0]]).GetSubstructMatches(Chem.MolFromSmarts('I'))
    mw = Chem.RWMol(Chem.MolFromSmiles(monomers[sequence[0]]))
    mw.ReplaceAtom(hits[0][0],Chem.Atom(17))
    Chem.SanitizeMol(mw)
    mw.CommitBatchEdit()
    polymer = Chem.AddHs(mw)
    info = Chem.AtomPDBResidueInfo()
    info.SetResidueName(str(chain_num) + sequence[0] + str(1))
    info.SetResidueNumber(1)
    [atom.SetMonomerInfo(info)  for  atom  in  polymer.GetAtoms()]
    Chem.SanitizeMol(polymer)
    
    for i in range(len(sequence))[1:]:
        if sequence[i] == 'A':
            A = Chem.MolFromSmiles(monomers['A'])
            A = Chem.AddHs(A)
            info = Chem.AtomPDBResidueInfo()
            info.SetResidueName(str(chain_num) + 'A' + str(i+1))
            info.SetResidueNumber(i+1)
            [atom.SetMonomerInfo(info)  for  atom  in  A.GetAtoms()]
            polymer = reaction.RunReactants((polymer, A))[0][0]
            Chem.SanitizeMol(polymer)
                
        elif sequence[i] == 'B':
            B = Chem.MolFromSmiles(monomers['B'])
            B = Chem.AddHs(B)
            info = Chem.AtomPDBResidueInfo()
            info.SetResidueName(str(chain_num) + 'B' + str(i+1))
            info.SetResidueNumber(i+1)
            [atom.SetMonomerInfo(info)  for  atom  in  B.GetAtoms()]
            polymer = reaction.RunReactants((polymer, B))[0][0]
            Chem.SanitizeMol(polymer)
        
        elif sequence[i] == 'S':
            S = Chem.MolFromSmiles(monomers['S'])
            S = Chem.AddHs(S)
            info = Chem.AtomPDBResidueInfo()
            info.SetResidueName(str(chain_num) + 'S' + str(i+1))
            info.SetResidueNumber(i+1)
            [atom.SetMonomerInfo(info)  for  atom  in  S.GetAtoms()]
            polymer = reaction.RunReactants((polymer, S))[0][0]
            Chem.SanitizeMol(polymer)
    
    if terminal == 'hydroxyl':
        hydrogen = Chem.MolFromSmiles('[H]')
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(str(chain_num) + sequence[0] + str(1))
        info.SetResidueNumber(1)
        [atom.SetMonomerInfo(info)  for  atom  in  hydrogen.GetAtoms()]
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('Cl'), hydrogen)[0]
        Chem.AddHs(polymer)
    elif terminal == 'carboxyl':
        carboxyl = Chem.MolFromSmiles('C(=O)[OH]')
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(str(chain_num) + sequence[0] + str(1))
        info.SetResidueNumber(1)
        [atom.SetMonomerInfo(info)  for  atom  in  carboxyl.GetAtoms()]
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('Cl'), carboxyl)[0]
    elif terminal == 'ester':
        carbon = Chem.MolFromSmiles('[CH3]')
        carbon = Chem.AddHs(carbon)
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(str(chain_num) + sequence[0] + str(1))
        info.SetResidueNumber(1)
        [atom.SetMonomerInfo(info)  for  atom  in  carbon.GetAtoms()]
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('Cl'), carbon)[0]
        #Chem.AddHs(polymer)
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('Cl'), Chem.MolFromSmiles('C'))[0]
    else:
        hydrogen = Chem.MolFromSmiles('[H]')
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(str(chain_num) + sequence[0] + str(1))
        info.SetResidueNumber(1)
        [atom.SetMonomerInfo(info)  for  atom  in  hydrogen.GetAtoms()]
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('Cl'), hydrogen)[0] #remove any excess iodine
    hydrogen = Chem.MolFromSmiles('[H]')
    info = Chem.AtomPDBResidueInfo()
    info.SetResidueName(str(chain_num) + sequence[-1] + str(len(sequence)))
    info.SetResidueNumber(len(sequence))
    [atom.SetMonomerInfo(info)  for  atom  in  hydrogen.GetAtoms()]
    polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('I'), hydrogen)[0] #remove any excess iodine
    Chem.SanitizeMol(polymer)
    return polymer



def build_linear_copolymer(sequence, 
                           monomer_a_smiles, 
                           monomer_b_smiles,
                           reaction=AllChem.ReactionFromSmarts('[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]')):
    """
    Constructs a linear co-polymer from the provided sequence of monomers.

    This function takes a sequence of monomers represented as 'A' and 'B', and the SMILES strings of two monomers.
    It constructs a co-polymer based on the sequence, using the provided reaction SMARTS for joining the monomers.
    The function returns the sanitized polymer and the percentage composition of each monomer in the polymer.

    Parameters
    ----------
    sequence : str
        A string representing the sequence of monomers. 'A' represents monomer_a and 'B' represents monomer_b.
    monomer_a_smiles : str
        The SMILES string of monomer A.
    monomer_b_smiles : str
        The SMILES string of monomer B.
    reaction : rdkit.Chem.rdChemReactions.ChemicalReaction, optional
        The reaction SMARTS used for joining the monomers. Defaults to '[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]',
        representing a condensation polymerisation.

    Returns
    -------
    tuple
        A tuple containing the following elements:
        - sanitized_polymer (rdkit.Chem.rdchem.Mol): The constructed and sanitized polymer.
        - percentage_monomer_a (float): The percentage composition of monomer A in the polymer.
        - percentage_monomer_b (float): The percentage composition of monomer B in the polymer.
    """
    # Initialize the polymer with an iodine blocker
    polymer = Chem.MolFromSmiles('OC(=O)I')
    A_count=0
    B_count=0
    A = Chem.MolFromSmiles(monomer_a_smiles)
    B = Chem.MolFromSmiles(monomer_b_smiles)
    # Build the polymer based on the sequence
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            polymer = reaction.RunReactants((polymer, A))[0][0]
            Chem.SanitizeMol(polymer)
            A_count+=1
        
        elif sequence[i] == 'B':
            polymer = reaction.RunReactants((polymer, B))[0][0]
            Chem.SanitizeMol(polymer)
            B_count+=1
    # Calculate the percentage composition of each monomer
    A_ratio = round((A_count/len(sequence))*100,2)
    B_ratio = round((B_count/len(sequence))*100,2)
    # Remove the iodine blocker
    polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('OC(=O)I'), Chem.MolFromSmarts('[O]'))[0]
    Chem.SanitizeMol(polymer)
    return polymer, A_ratio, B_ratio



def PDI(chains):
    """
    Calculates the Polydispersity Index (PDI), number-average molecular weight (Mn),
    and weight-average molecular weight (Mw) of a list of chains.

    Parameters
    ----------
    chains : list
        A list of molecular chains. Each chain is represented as an RDKit molecule object.

    Returns
    -------
    tuple
        A tuple containing the following elements:
        - PDI (float): The Polydispersity Index, which is the ratio of Mw to Mn.
        - Mn (float): The number-average molecular weight.
        - Mw (float): The weight-average molecular weight.
    """
    # Calculate the molecular weights of the chains
    mw_list = [ExactMolWt(chain) for chain in chains]
    # Count instances of each molecular weight
    mw_counts = Counter(mw_list)
    # Calculate Mn 
    total_chains = sum(mw_counts.values())  # Total number of chains 
    Mn = sum(mw * count for mw, count in mw_counts.items()) / total_chains
    # Calculate Mw
    Mw = sum(mw**2 * count for mw, count in mw_counts.items()) / sum(mw * count for mw, count in mw_counts.items())
    # Calculate PDI - Polydispersity Index
    PDI = Mw / Mn
    
    return PDI, Mn, Mw
    



def blockiness_gen(sequence, wrt='A'):
    """
    Calculate the blockiness and average block length of a co-polymer sequence.

    This function calculates the blockiness of a co-polymer sequence by counting the occurrences of 'BB' and 'BA' or 'AB' in the sequence.
    It also calculates the average block length of 'A' and 'B' monomers in the sequence.

    Parameters
    ----------
    sequence : str
        A string representing the co-polymer sequence. 'A' represents one type of monomer and 'B' represents another type.

    Returns
    -------
    tuple
        A tuple containing the following elements:
        - blockiness (float): The blockiness of the co-polymer sequence. Calculated as the ratio of 'BB' to 'BA' or 'AB'.
        - block_length_A (float): The average block length of 'A' in the sequence.
        - block_length_B (float): The average block length of 'B' in the sequence.

    Notes
    -----
    If the sequence does not contain both 'A' and 'B', the function returns a string indicating that the molecule is not a co-polymer.
    """
    if wrt == 'A':
    
        if 'A' in sequence and 'B' in sequence: #Check if sequence is a co-polymer
            AB = sequence.count('AB')
            BB = sequence.count('BB')
            BA = sequence.count('BA')
            AA = sequence.count('AA')
            if 'BA' in sequence:
                blockiness = AA/BA
            else:
                blockiness = AA/AB
            #Calculate block length B
            block_list_B = [x for x in sequence.split('A') if x!='']
            block_length_B = mean([len(b) for b in block_list_B])
            #Calculate block length A
            block_list_A = [x for x in sequence.split('B') if x!='']
            block_length_A = mean([len(b) for b in block_list_A])
            return blockiness, block_length_B, block_length_A

        else:
            return 'Molecule is not a co-polymer, no blockiness calculation performed', 0, len(sequence)
    
    elif wrt == 'B':
       
        if 'A' in sequence and 'B' in sequence:
            AB = sequence.count('AB')
            BB = sequence.count('BB')
            BA = sequence.count('BA')
            AA = sequence.count('AA')
            if 'AB' in sequence:
                blockiness = BB/AB
            else:
                blockiness = BB/BA
            #Calculate block length A
            block_list_A = [x for x in sequence.split('B') if x!='']
            block_length_A = mean([len(b) for b in block_list_A])
            #Calculate block length B
            block_list_B = [x for x in sequence.split('A') if x!='']
            block_length_B = mean([len(b) for b in block_list_B])
            return blockiness, block_length_A, block_length_B

        else:  
            return 'Molecule is not a co-polymer, no blockiness calculation performed', 0, len(sequence)
    

    else:
        raise ValueError("wrt parameter must be 'A' or 'B'")

def calculate_box_components(chains, monomers, sequence, salt_concentration = 0.0 * unit.mole / unit.liter, residual_monomer = 0.00, solvated=False):
    """
    Calculates the components required to construct a simulation box for a given set of molecular chains.
    
    This function determines the quantity of each molecule type required, considering the salt concentration
    and residual monomer concentration. It is adapted from the OpenFF Toolkit Packmol wrapper's solvate_topology function.
    
    Parameters:
    -----------
    chains : list
        A list of molecular chains to be included in the simulation box.
    sequence : str
        A string representing the sequence of the molecular chains. 'G' and 'L' represent different types of monomers.
    salt_concentration : float, optional
        The desired salt concentration in the simulation box. Defaults to 0 M.
    residual_monomer : float, optional
        The desired residual monomer concentration in the simulation box. Defaults to 0.00%.
    solvated : bool, optional
        Indicates whether the system contains water. Defaults to False.
    
    Returns:
    --------
    tuple
        A tuple containing the following elements:
        - molecules (list): A list of molecules to be included in the simulation box.
        - number_of_copies (list): A list indicating the quantity of each molecule to be included in the simulation box.
        - topology (openff.toolkit.topology.Topology): The topology of the simulation box.
        - box_vectors (numpy.ndarray): The vectors defining the dimensions of the simulation box.
    
    Notes:
    ------
    This function is adapted from the OpenFF Toolkit Packmol wrapper's solvate_topology function.
    """
    from openff.toolkit.topology import Molecule, Topology
    from openff.interchange.components._packmol import UNIT_CUBE, pack_box, RHOMBIC_DODECAHEDRON, solvate_topology
    from openff.interchange.components._packmol import _max_dist_between_points, _compute_brick_from_box_vectors, _center_topology_at
    import warnings

    #Create molecules for the purpose of mass calculation
    #Water
    water = Molecule.from_smiles('O')
    water.generate_unique_atom_names()
    water.generate_conformers()
    water_mass = sum([atom.mass for atom in water.atoms])
    
    #Sodium
    na = Molecule.from_smiles('[Na+]')
    na.generate_unique_atom_names()
    na.generate_conformers()
    
    #Chloride
    cl = Molecule.from_smiles('[Cl-]')
    cl.generate_unique_atom_names()
    cl.generate_conformers()
    
    nacl_mass = sum([atom.mass for atom in na.atoms]) + sum(
    [atom.mass for atom in cl.atoms],)
    # Create a topology from the chains
    topology = Topology.from_molecules(chains)
    nacl_conc=salt_concentration
    padding= 0.1 * unit.nanometer
    box_shape= UNIT_CUBE
    target_density= 1.0 * unit.gram / unit.milliliter

    
    # Compute box vectors from the solute length and requested padding
    if chains[0].n_conformers == 0:
        raise ValueError("The solvate_topology function requires that the solute has at least one conformer.")
    solute_length = max(_max_dist_between_points(chains[i].to_topology().get_positions()) for i in range(len(chains)))
    image_distance = solute_length + padding * 2
    box_vectors = box_shape * image_distance
    
    # Compute target masses of solvent
    box_volume = np.linalg.det(box_vectors.m) * box_vectors.u**3
    target_mass = box_volume * target_density
    solvent_mass = target_mass - sum(sum([atom.mass for atom in molecule.atoms]) for molecule in topology.molecules)
    
    # Compute the number of NaCl to add from the mass and concentration
    nacl_mass_fraction = (nacl_conc * nacl_mass) / (55.5 * unit.mole / unit.liter * water_mass)
    nacl_to_add = ((solvent_mass * nacl_mass_fraction) / nacl_mass).m_as(unit.dimensionless).round()
    if solvated:
        water_to_add = int(round((solvent_mass - nacl_mass) / water_mass).m_as(unit.dimensionless).round())
    else:
        water_to_add = 0
    
    # Neutralise the system by adding and removing salt
    solute_charge = sum([molecule.total_charge for molecule in topology.molecules])
    na_to_add = int(round(np.ceil(nacl_to_add - solute_charge.m / 2.0)))
    cl_to_add = int(round(np.floor(nacl_to_add + solute_charge.m / 2.0)))
    
    rolling_mass=0
    for m in topology.molecules:
        rolling_mass += sum(atom.mass for atom in m.atoms)
    rolling_mass += nacl_mass * nacl_to_add

    # Calculate residual monomer to add
    mass_to_add = (rolling_mass.magnitude/100-residual_monomer) * residual_monomer
    if mass_to_add < 90:
        warnings.warn('Residual monomer calculation may not be accurate for small systems, as the residual monomer mass can be lower than the weight of a single monomer. Please check the output by querying the value of residual_monomer_actual')
    
    if 'A' in sequence and 'B' in sequence:
        A_rd = Chem.MolFromSmiles(monomers[0])
        A_rd = Chem.AddHs(A_rd)
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName('A' + str(1))
        info.SetResidueNumber(1)
        [atom.SetMonomerInfo(info)  for  atom  in  A_rd.GetAtoms()]
        A = Molecule.from_rdkit(A_rd)
        A.name = 'A_residual'
        A_mass = sum([atom.mass for atom in A.atoms])
        B_rd = Chem.MolFromSmiles(monomers[1])
        B_rd = Chem.AddHs(B_rd)
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName('B' + str(1))
        info.SetResidueNumber(1)
        [atom.SetMonomerInfo(info)  for  atom  in  B_rd.GetAtoms()]
        B = Molecule.from_rdkit(B_rd)
        B_mass = sum([atom.mass for atom in B.atoms])
        B.name = 'B_residual'
        for r in range(0,100):
            if (r * A_mass.magnitude) + (r * B_mass.magnitude) <= mass_to_add:
                A_to_add = r
                B_to_add = r
            else:
                break
        residual_monomer_actual = ((A_to_add * A_mass.magnitude + B_to_add * B_mass.magnitude) / rolling_mass.magnitude) *100
        molecules = [water, na, cl, A, B]
        number_of_copies=[water_to_add, na_to_add, cl_to_add, A_to_add, B_to_add]
    
    elif 'A' in sequence and 'B' not in sequence:
        A_rd = Chem.MolFromSmiles(monomers[0])
        A_rd = Chem.AddHs(A_rd)
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName('A' + str(1))
        info.SetResidueNumber(1)
        [atom.SetMonomerInfo(info)  for  atom  in  A_rd.GetAtoms()]
        A = Molecule.from_rdkit(A_rd)
        A_mass = sum([atom.mass for atom in A.atoms])
        A.name = 'A_residual'
        B = Molecule.from_smiles('C')
        B.name = 'B_residual'
        for r in range(0,100):
            if r * A_mass.magnitude <= mass_to_add:
                A_to_add = r
            else:
                break
        B_to_add = 0
        residual_monomer_actual = ((A_to_add * A_mass.magnitude) / rolling_mass.magnitude) * 100
        molecules = [water, na, cl, A, B]
        number_of_copies=[water_to_add, na_to_add, cl_to_add, A_to_add, B_to_add]
    return molecules, number_of_copies, topology, box_vectors, residual_monomer_actual

def introduce_stereoisomers(stereo_monomer, instance, seq):
    """
    Introduce stereoisomers into a polymer sequence by replacing a specified percentage of 'A' monomers with 'S'.

    This function replaces a specified percentage of 'A' monomers with 'S' in the given sequence. The replacements
    are made in pairs of 'A' monomers to ensure stereoisomerism.

    Parameters
    ----------
    stereo_monomer : str
        The monomer to be replaced with its stereoisomer (e.g., 'A').
    instance : float
        The fraction of 'stereo_monomer' monomers to be replaced with stereoisomer (e.g., 0.5 for 50%).
    seq : str
        The original polymer sequence.

    Returns
    -------
    str
        The modified sequence with stereoisomers introduced.

    """
    assert stereo_monomer in seq, f"Monomer {stereo_monomer} not found in sequence"
    assert instance <= 1 and instance > 0, "Instance must be a fraction between 0 and 1"
    assert instance > 0, "Instance must be greater than 0"
    assert type(stereo_monomer) == str, "Monomer must be a string"
    num_replacements = round(seq.count(stereo_monomer) * instance)
    seq_list = list(seq)

    # Initialize a counter for replacements
    replacement_count = 0
    i = 0
    while i < len(seq_list) - 1:
        if seq_list[i] == 'A' and seq_list[i + 1] == 'A':
            seq_list[i + 1] = 'S'
            replacement_count += 1
            if replacement_count >= num_replacements:
                break
            i += 1
        i += 1
    
    modified_seq = ''.join(seq_list)
    return modified_seq

#Class object for generic polymer system

class polymer_system:
    try:
        from openeye import oechem
    except:
        import warnings
        warnings.warn("OpenEye is not installed. You will not be able to use OpenEye Toolkits for conformer generation.")
    from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
    from functools import reduce
    from statistics import mean
    from rdkit.Chem.Descriptors import ExactMolWt
    from openff.interchange import Interchange
    from openff.interchange.components._packmol import UNIT_CUBE, pack_box
    from swiftpol.build import build_polymer, PDI, blockiness_gen, calculate_box_components, introduce_stereoisomers
    from openff.units import unit
    from rdkit.Chem import AllChem

    def __init__(self, 
                 monomer_list, 
                 reaction, 
                 length_target, 
                 num_chains, 
                 stereoisomerism_input=None,
                 terminals='standard', 
                 perc_A_target=100, 
                 blockiness_target=[1.0, 'A'], 
                 copolymer=False, 
                 acceptance = 10
                 ):
        """
        Initialize the polymer system and build the polymer chains.

        **Parameters:**
        ------------

        monomer_list (list): List of monomers to be used in the polymerization.

        reaction (str): The type of reaction to be used for polymerization.

        length_target (float): The target length of the polymer chains.

        num_chains (int): The number of polymer chains to be generated.

        stereoisomerism_input (tuple, optional): A tuple containing the monomer, instance fraction (e.g. 0.5 for 50% stereoisomer), and SMILES string of the stereoisomer to be introduced. Default is None.

        terminals (str, optional): The type of terminal groups to be used. Default is 'standard'.

        perc_A_target (float, optional): The target percentage of monomer A in the copolymer. Default is 100.

        blockiness_target (float, optional): The target blockiness of the copolymer, and indication of method calculation. Default is 1.0, with reference to 'A' monomer linkages.

        copolymer (bool, optional): Flag to indicate if the system is a copolymer. Default is False.

        acceptance = % deviation of blockiness and A percentage from target values. Default is 10%


        **Attributes:**
        ---------------

        length_target (float): The target length of the polymer chains.

        terminals (str): The type of terminal groups used.

        blockiness_target (float): The target blockiness of the copolymer.

        A_target (float): The target percentage of monomer A in the copolymer.

        chains (list): List of polymer chains as OpenFF Molecule objects.

        chain_rdkit (list): List of polymer chains as RDKit molecule objects.

        lengths (list): List of lengths of the polymer chains.

        perc_A_actual (list): List of actual percentages of monomer A in the polymer chains.

        B_block_length (float): The average block length of monomer B in the copolymer.

        A_block_length (float): The average block length of monomer A in the copolymer.

        blockiness_list (list): List of blockiness values for the polymer chains.

        mean_blockiness (float): The mean blockiness of the polymer chains.

        mol_weight_average (float): The average molecular weight of the polymer chains.

        PDI (float): The polydispersity index of the polymer chains.

        Mn (float): The number-average molecular weight of the polymer chains.

        Mw (float): The weight-average molecular weight of the polymer chains.

        num_chains (int): The number of polymer chains generated.

        length_average (float): The average length of the polymer chains.

        min_length (float): The minimum length of the polymer chains.

        max_length (float): The maximum length of the polymer chains.

        """
        self.length_target = length_target
        self.terminals = terminals
        perc_A_actual = []
        if copolymer==True:
            self.blockiness_target = blockiness_target[0]
            self.A_target = perc_A_target
            
            def spec(sequence): #Define limits of A percentage and blockiness from input
                acceptance_dec = acceptance/100
                actual_A = (sequence.count('A')/len(sequence))*100
                blockiness = blockiness_gen(sequence, blockiness_target[1])[0]
                return actual_A > perc_A_target*(1-acceptance_dec) and actual_A < perc_A_target*(1+acceptance_dec) and blockiness>blockiness_target[0]*(1-acceptance_dec) and blockiness<blockiness_target[0]*(1+acceptance_dec)
            
            blockiness_list = []
            out_of_spec = 0
            BBL = []
            ABL = []
        chains = []
        chains_rdkit = []
        lengths = []
        
        self.monomers = [mono.replace("[I]", "") for mono in monomer_list]
        self.reaction = reaction
        reaction = AllChem.ReactionFromSmarts(reaction)
        if stereoisomerism_input is not None:
            stereo_monomer, instance, new_smiles = stereoisomerism_input
            monomer_list.append(new_smiles)
        #First round of building - copolymer
        if copolymer==True:
            for n in range(num_chains):
                length_actual = np.random.normal(length_target, 0.5)
                sequence = reduce(lambda x, y: x + y, np.random.choice(['A', 'B'], size=(int(length_actual),), p=[perc_A_target/100,1-(perc_A_target/100)]))
                blockiness = blockiness_gen(sequence, blockiness_target[1])[0]
                if spec(sequence)==True:
                    if stereoisomerism_input is not None:
                        sequence_stereo = introduce_stereoisomers(stereo_monomer, instance, sequence)
                        pol = build_polymer(sequence=sequence_stereo, monomer_list = monomer_list, reaction = reaction, terminal=terminals, chain_num=n+1)
                    else:
                        pol = build_polymer(sequence=sequence, monomer_list = monomer_list, reaction = reaction, terminal=terminals, chain_num=n+1)
                    lengths.append(int(length_actual))
                    chains_rdkit.append(pol)
                    chain = Molecule.from_rdkit(pol)
                    chains.append(chain)
                    perc_A_actual.append((sequence.count('A')/len(sequence))*100)
                    blockiness_list.append(blockiness)
                    BBL.append(blockiness_gen(sequence, blockiness_target[1])[1])
                    ABL.append(blockiness_gen(sequence, blockiness_target[1])[2])
                else:
                    out_of_spec +=1
                #Second round of building
                while out_of_spec >0:
                    length_actual = np.random.normal(length_target, 0.5)
                    sequence = reduce(lambda x, y: x + y, np.random.choice(['A', 'B'], size=(int(length_actual),), p=[perc_A_target/100,1-(perc_A_target/100)]))
                    blockiness = blockiness_gen(sequence, blockiness_target[1])[0]
                    if spec(sequence)==True:
                        if stereoisomerism_input is not None:
                            sequence_stereo = introduce_stereoisomers(stereo_monomer, instance, sequence) 
                            pol = build_polymer(sequence=sequence_stereo, monomer_list = monomer_list, reaction = reaction, terminal=terminals, chain_num=n+1)
                        else:
                            pol = build_polymer(sequence=sequence, monomer_list = monomer_list, reaction = reaction, terminal=terminals, chain_num=n+1)                   
                        lengths.append(int(length_actual))
                        chains_rdkit.append(pol)
                        chain = Molecule.from_rdkit(pol)
                        chains.append(chain)
                        perc_A_actual.append((sequence.count('A')/len(sequence))*100)
                        blockiness_list.append(blockiness)
                        BBL.append(blockiness_gen(sequence, blockiness_target[1])[1])
                        ABL.append(blockiness_gen(sequence, blockiness_target[1])[2])
                        out_of_spec-=1

                self.B_block_length = mean(BBL)
                self.A_block_length = mean(ABL)
                self.blockiness_list = blockiness_list
                self.mean_blockiness = mean(blockiness_list)
                self.perc_A_actual = perc_A_actual
                self.A_actual = mean(perc_A_actual)
        else:
            for n in range(num_chains):
                length_actual = np.random.normal(length_target, 0.5)
                sequence = reduce(lambda x, y: x + y, np.random.choice(['A', 'B'], size=(int(length_actual),), p=[perc_A_target/100,1-(perc_A_target/100)]))
                if stereoisomerism_input is not None:
                    sequence_stereo = introduce_stereoisomers(stereo_monomer, instance, sequence) 
                    pol = build_polymer(sequence=sequence_stereo, monomer_list = monomer_list, reaction = reaction, terminal=terminals, chain_num=n+1)
                else:
                    pol = build_polymer(sequence=sequence, monomer_list = monomer_list, reaction = reaction, terminal=terminals, chain_num=n+1)                   
                lengths.append(int(length_actual))
                chains_rdkit.append(pol)
                chain = Molecule.from_rdkit(pol)
                chains.append(chain)
                perc_A_actual.append((sequence.count('A')/len(sequence))*100)
            self.B_block_length = None
            self.A_block_length = None
            self.blockiness_list = None
            self.mean_blockiness = None
            self.perc_A_actual = None
            self.A_actual = None
        self.sequence = sequence
        self.chains = chains
        for i in range(len(self.chains)):
            self.chains[i].name = 'chain' + str(i+1)
        self.chain_rdkit = chains_rdkit
        self.mol_weight_average = round(mean([ExactMolWt(c) for c in chains_rdkit]),2)
        self.PDI, self.Mn, self.Mw = PDI(chains_rdkit)
        self.num_chains = len(chains)
        self.A_actual = mean(perc_A_actual)
        self.perc_A_actual = perc_A_actual
        self.length_average = mean(lengths)
        self.lengths = lengths
        self.min_length = min(lengths)
        self.max_length = max(lengths)
        print('System built!, size =', self.num_chains)
        
    def __repr__(self):
                
        description = (f"SwiftPol ensemble of size {self.num_chains}, "
                       f"average chain length = {self.length_average}-mers, PDI = {self.PDI}")
        return description

    def generate_conformers(self):
        """
        Generate conformers for each polymer chain in the system.

        This method uses the OpenFF toolkit OpenEye Wrapper to generate conformers for each polymer chain in the system.
        It first checks if the OpenEye toolkit is licensed and available. If it is, it uses the OpenEyeToolkitWrapper
        to generate conformers. Otherwise, it falls back to using the RDKitToolkitWrapper. Each chain is processed
        to generate a single conformer, and unique atom names are assigned to each chain.

        Raises
        ------
        ImportError
            If neither RDKit nor OpenEye toolkits are available.
        """

        from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
        # Generate conformers using OpenFF toolkit wrapper
        for chain in self.chains:
            num = self.chains.index(chain)
            object = RDKitToolkitWrapper()
            if oechem_imported:
                if oechem.OEChemIsLicensed():
                    object = OpenEyeToolkitWrapper()
            object.generate_conformers(molecule=chain, n_conformers=1)
            chain.generate_unique_atom_names()
            self.chains[num] = chain
    
    def charge_system(self, charge_scheme):
        """
        Assign partial charges to each polymer chain in the system.

        This method uses one of AM1-BCC, Espaloma, or OpenFF NAGL to assign partial charges to each polymer chain in the system.
        It iterates over each chain in the `self.chains` list and assigns partial charges to the chain.

        Parameters
        ----------
        charge_scheme : str
            The charge assignment scheme to use. Options are 'AM1_BCC', 'espaloma', or 'NAGL'.

        Raises
        ------
        ImportError
            If the selected toolkit is not available.
        """
        from swiftpol.parameterize import charge_openff_polymer
        for chain in self.chains:
            chain.partial_charges = charge_openff_polymer(chain, charge_scheme)
    def export_to_csv(self, filename):
        """
        Export all the instances in polymer_system.__init__() into a pandas DataFrame and save it as a CSV file.

        Parameters
        ----------
        filename : str
            The name of the CSV file to save the data.
        """
        data = {
            'monomers': [self.monomers],
            'length_target': [self.length_target],
            'terminals': [self.terminals],
            'num_chains': [self.num_chains],
            'mol_weight_average': [self.mol_weight_average],
            'PDI': [self.PDI],
            'Mn': [self.Mn],
            'Mw': [self.Mw],
            'sequence': [self.sequence],
            'B_block_length': [self.B_block_length],
            'A_block_length': [self.A_block_length],
            'blockiness_list': [self.blockiness_list],
            'mean_blockiness': [self.mean_blockiness],
            'perc_A_actual': [self.perc_A_actual],
            'A_actual': [self.A_actual]
        }

        df = pd.DataFrame(data)
        df.to_csv(filename, index=False)

    def pack_solvated_system(self, salt_concentration=0.0 * unit.mole / unit.liter, residual_monomer=0.00):
        """
        Pack a solvated system using Packmol functions, and the OpenFF Packmol wrapper.
    
        This method uses Packmol to build a solvated system by packing molecules into a simulation box.
        It considers the salt concentration and residual monomer concentration to determine the quantity of each molecule type required.
    
        Parameters
        ----------
        salt_concentration : openff.units.Quantity, optional
            The desired salt concentration in the simulation box. Default is 0.0 M.
        residual_monomer : float, optional
            The desired residual monomer concentration in the simulation box. Default is 0.00.
    
        Returns
        -------
        openff.toolkit.topology.topology.Topology
            An Interchange object representing the packed solvated system.

        Notes
        -----
        This function uses the OpenFF Interchange Packmol wrapper to pack the molecules into the simulation box.
        It removes any molecules with zero copies before packing, to avoid packmol errors.
        Assigns % residual monomer value to ensemble under self.residual_monomer_actual

        """
        from openff.interchange.components._packmol import pack_box, UNIT_CUBE
        from swiftpol.build import calculate_box_components
    
        molecules, number_of_copies, topology, box_vectors, residual_monomer_actual = calculate_box_components(
            self.chains, self.monomers, self.sequence, salt_concentration, residual_monomer
        )
        molecules = [molecules[i] for i in range(len(number_of_copies)) if number_of_copies[i] != 0]
        number_of_copies = [num for num in number_of_copies if num != 0]
        self.residual_monomer_actual = residual_monomer_actual
        if topology.n_molecules == 1:
            return pack_box(molecules=molecules,
                            number_of_copies=number_of_copies,
                            solute=topology,
                            box_vectors=box_vectors,
                            box_shape=UNIT_CUBE,
                            center_solute = 'BRICK')
        else:
            return pack_box(molecules = molecules + self.chains,
                            number_of_copies = number_of_copies+[1 for i in range(len(self.chains))],
                            box_vectors = box_vectors,
                            tolerance = 1*unit.angstrom)

    def generate_polyply_files(self, residual_monomer=0.00, residual_oligomer=0.00):
        """
        Generate input files for Polyply from the system.

        This method generates the input files required for Polyply (https://github.com/marrink-lab/polyply_1.0/) from the system.

        Parameters
        ----------
        residual_monomer : float, optional
            The desired residual monomer concentration in the simulation box. Default is 0.00.
        residual_oligomer : float, optional
            The desired residual oligomer concentration in the simulation box. Default is 0.00.

        Returns
        -------
        tuple
            A tuple containing the paths to the .gro, .top, and .pdb files generated for Polyply and GROMACS insert-molecules input.

        Notes
        -----
        This function uses OpenFF Interchange to generate the input files for Polyply.

        Raises
        ------
        UserWarning
            If partial charges are not assigned to the system, processing large systems may raise errors from OpenFF-Interchange.
        """
        from swiftpol.build import calculate_box_components
        from openff.interchange import Interchange
        from openff.toolkit import ForceField
        import warnings
        box_vectors = calculate_box_components(self.chains, self.monomers, self.sequence, 0.0 * unit.mole / unit.liter, 0.0)[3]
        molecules, number_of_copies, residual_monomer_actual, residual_oligomer_actual = self.calculate_residuals(residual_monomer, residual_oligomer)
        mol_pdb_files_dest = []
        for i in molecules:
            string_i = str(molecules.index(i)) + '.pdb'
            mol_pdb_files_dest.append(string_i)
            i.generate_conformers(n_conformers=1)
            i.to_file(string_i, file_format='pdb')
        self.residual_monomer_actual = residual_monomer_actual
        topology = Topology.from_molecules(self.chains+[molecules[i] for i in range(len(number_of_copies)) if number_of_copies[i] != 0])
        if self.chains[0].partial_charges is None:
            warnings.warn('Partial charges may not be assigned to the system. Processing large systems may raise errors from OpenFF-Interchange', UserWarning)
            interchange = Interchange.from_smirnoff(
                topology=topology,
                force_field=ForceField("openff-2.2.0.offxml"),
                box=box_vectors
            )
        else:
            interchange = Interchange.from_smirnoff(
                topology=topology,
                force_field=ForceField("openff-2.2.0.offxml"),
                charge_from_molecules=[i for i in self.chains],
                box=box_vectors
            )
        interchange.to_gromacs('swiftpol_output')
        return_tuple = ('swiftpol_output.gro', 'swiftpol_output.top')
        for pdb_file in mol_pdb_files_dest:
            return_tuple += (pdb_file,)
        print(f'Polyply input files generated! Saved at {return_tuple}')
        return return_tuple

    def calculate_residuals(self, residual_monomer = 0, residual_oligomer = 0):
        """
        Generate residual monomer and oligomer molecules, and molecule counts.

        Parameters
        ----------
        residual_monomer : float
            The desired residual monomer concentration in the simulation box. Default is 0.00.
        residual_oligomer : float
            The desired residual oligomer concentration in the simulation box. Default is 0.00.
            
        Returns
        -------
        tuple
        A tuple containing the following elements:
            - A list of OpenFF Molecule objects representing the residual monomer and oligomer molecules.
            - A list of molecule counts, corresponding with the molecule list
            - Actual value for % residual monomer
            - Actual value for % residual oligomer

        Notes
        -----
        EXPERIMENTAL CAPABILITY. Proceed with caution

        This function works independently of, and is unrelated to, the build.calculate_box_components function.
        
        Raises
        ------
        UserWarning - If the residual monomer concentration is close to or lower than the weight of a single monomer.

        
        """
        from functools import reduce
        from swiftpol import build
        import numpy as np
        from openff.toolkit import Topology, Molecule
        from rdkit import Chem
        import warnings
        topology = Topology.from_molecules(self.chains)
        rolling_mass=0
        monomers = self.monomers
        sequence = self.sequence
        for m in topology.molecules:
            rolling_mass += sum(atom.mass for atom in m.atoms)
        if rolling_mass.magnitude < 1000:
            warnings.warn('Residual monomer/oligomer calculation may not be accurate for small systems, as the residual mass can be close to or lower than the weight of a single monomer. Please check the output by querying the value of residual_monomer_actual')
        number_of_copies = []
        # residual to add
        monomer_to_add = (rolling_mass.magnitude/100 - residual_monomer - residual_oligomer) * residual_monomer
        oligo_to_add = (rolling_mass.magnitude/100 - residual_monomer - residual_oligomer) * residual_oligomer
        if 'A' in sequence and 'B' in sequence:
            A_rd = Chem.MolFromSmiles(monomers[0])
            A_rd = Chem.AddHs(A_rd)
            info = Chem.AtomPDBResidueInfo()
            info.SetResidueName('A' + str(1))
            info.SetResidueNumber(1)
            [atom.SetMonomerInfo(info)  for  atom  in  A_rd.GetAtoms()]
            A = Molecule.from_rdkit(A_rd)
            A.name = 'A_residual'
            A_mass = sum([atom.mass for atom in A.atoms])
            B_rd = Chem.MolFromSmiles(monomers[1])
            B_rd = Chem.AddHs(B_rd)
            info = Chem.AtomPDBResidueInfo()
            info.SetResidueName('B' + str(1))
            info.SetResidueNumber(1)
            [atom.SetMonomerInfo(info)  for  atom  in  B_rd.GetAtoms()]
            B = Molecule.from_rdkit(B_rd)
            B_mass = sum([atom.mass for atom in B.atoms])
            B.name = 'B_residual'
            for r in range(1,100):
                if (r * A_mass.magnitude) + (r * B_mass.magnitude) >= monomer_to_add:
                    mon_count = r
                    break
            A_to_add = mon_count
            B_to_add = mon_count

        elif 'A' in sequence and 'B' not in sequence:
            A_rd = Chem.MolFromSmiles(monomers[0])
            A_rd = Chem.AddHs(A_rd)
            info = Chem.AtomPDBResidueInfo()
            info.SetResidueName('A' + str(1))
            info.SetResidueNumber(1)
            [atom.SetMonomerInfo(info)  for  atom  in  A_rd.GetAtoms()]
            A = Molecule.from_rdkit(A_rd)
            A_mass = sum([atom.mass for atom in A.atoms])
            A.name = 'A_residual'
            B = Molecule.from_smiles('C')
            B.name = 'B_residual'
            for r in range(1,100):
                if (r * A_mass.magnitude) >= monomer_to_add:
                    mon_count = r
                    break

            A_to_add = mon_count
            B_to_add = 0
        
       
        oligomers = []
        if 'A' in sequence and 'B' in sequence:
            for i in range(1000):
                oligo_seq = reduce(lambda x, y: x + y, np.random.choice(['A', 'B'], size=(int(self.length_target * 0.1)), p=[self.A_target/100,1-(self.A_target/100)]))
                monomer_list = [mono+'[I]' for mono in self.monomers]
                oligomer_rd = build.build_polymer(oligo_seq, 
                            monomer_list=monomer_list, 
                            reaction=AllChem.ReactionFromSmarts(self.reaction))
                oligomer_rd = Chem.AddHs(oligomer_rd)
                info = Chem.AtomPDBResidueInfo()
                info.SetResidueName('O' + str(i+1))
                info.SetResidueNumber(1)
                [atom.SetMonomerInfo(info)  for  atom  in  oligomer_rd.GetAtoms()]
                oligomer = Molecule.from_rdkit(oligomer_rd)
                oligo_mass = 0
                oligomers_new = oligomers + [oligomer]
                for i in oligomers_new:
                    oligo_mass += sum(atom.mass for atom in i.atoms)
                if oligo_mass.magnitude <= oligo_to_add:
                    oligomers.append(oligomer)
                else:
                    break
            monomer_mass = (A_to_add * A_mass.magnitude) + (B_to_add * B_mass.magnitude)
            oligomer_mass = 0 * unit.dalton
            for i in oligomers:
                oligomer_mass += sum(atom.mass for atom in i.atoms)
            residual_monomer_actual = (monomer_mass / (rolling_mass.magnitude + monomer_mass + oligomer_mass.magnitude)) * 100
            residual_oligomer_actual = (oligomer_mass.magnitude / (rolling_mass.magnitude + monomer_mass + oligomer_mass.magnitude)) * 100
            
            molecules = [A, B] + oligomers
            number_of_copies = [A_to_add, B_to_add] + [1 for c in range(len(oligomers))]

        elif 'A' in sequence and 'B' not in sequence:
            for i in range(1000):
                oligo_seq = reduce(lambda x, y: x + y, np.random.choice(['A'], size=(int(self.length_target * 0.1)), p=[self.A_target/100,1-(self.A_target/100)]))
                monomer_list = [mono+'[I]' for mono in self.monomers]
                oligomer_rd = build.build_polymer(oligo_seq, 
                            monomer_list=monomer_list, 
                            reaction=AllChem.ReactionFromSmarts(self.reaction))
                oligomer_rd = Chem.AddHs(oligomer_rd)
                info = Chem.AtomPDBResidueInfo()
                info.SetResidueName('O' + str(i+1))
                info.SetResidueNumber(1)
                [atom.SetMonomerInfo(info)  for  atom  in  oligomer_rd.GetAtoms()]
                oligomer = Molecule.from_rdkit(oligomer_rd)
                oligomer.name = 'oligo' + str(i+1)
                oligo_mass = 0
                oligomers_new = oligomers + [oligomer]
                for i in oligomers_new:
                    oligo_mass += sum(atom.mass for atom in i.atoms)
                if oligo_mass.magnitude <= oligo_to_add:
                    oligomers.append(oligomer)
                else:
                    break

            monomer_mass = (A_to_add * A_mass.magnitude) + (B_to_add * B_mass.magnitude)
            oligomer_mass = 0 * unit.dalton
            for i in oligomers:
                oligomer_mass += sum(atom.mass for atom in i.atoms)
            residual_monomer_actual = (monomer_mass / (rolling_mass.magnitude + monomer_mass + oligomer_mass.magnitude)) * 100
            residual_oligomer_actual = (oligomer_mass.magnitude / (rolling_mass.magnitude + monomer_mass + oligomer_mass.magnitude)) * 100
            
            molecules = [A, B] + oligomers
            number_of_copies = [A_to_add, B_to_add] + [1 for c in range(len(oligomers))]
    
        return molecules, number_of_copies, residual_monomer_actual, residual_oligomer_actual