#Init
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import random
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
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
from swiftpol.build import PDI


from functools import reduce
from statistics import mean
from rdkit.Chem.Descriptors import ExactMolWt
from openff.interchange import Interchange
from swiftpol import build, demo


#Ring opening polymerisation
def build_PLGA_ring(sequence, 
                    reaction = AllChem.ReactionFromSmarts('[I:1][O:2].[I:3][C:4]>>[C:4][O:2].[I:3][I:1]'),
                    terminal='hydroxyl'):
    """
    Build a PLGA co-polymer of a specified diblock sequence and return the sanitized polymer.

    This function takes a list of up to 2 monomers to create a co-polymer. It uses the cyclic esters lactide and glycolide as constituent monomers.
    The function joins monomers using the specified reaction SMARTS and returns the constructed PLGA macromolecule along with the percentage composition of lactide and glycolide.

    Parameters
    ----------
    sequence : str
        A string representing the sequence of monomers. 'L' represents lactide and 'G' represents glycolic acid.
        The sequence must be assembled as blocks of 2 monomers (e.g., 'LLGGLLGG').
    monomer_list : list
        A list of RDKit.Chem.Mol objects representing the monomers. Each monomer must contain I-I bridging the point of esterification.
    reaction : rdkit.Chem.rdChemReactions.ChemicalReaction
        An RDKit reaction object specifying the joining of 2 iodinated compounds into an ester.

    Returns
    -------
    tuple
        A tuple containing the following elements:
        - plga_macromolecule (rdkit.Chem.rdchem.Mol): The constructed and sanitized PLGA macromolecule.
        - lactide_ratio (float): The percentage composition of lactide in the polymer.
        - glycolide_ratio (float): The percentage composition of glycolide in the polymer.
    """
    # Function implementation here
    ring_smiles = ['O1C(=O)C[I+][I+]OC(=O)C1', 'C[C@@H]1[I+][I+]OC(=O)[C@H](C)OC1=O', 'C[C@@H]1O[I+][I+]C(=O)[C@@H](C)OC1=O', 'C[C@H]1O[I+][I+]C(=O)[C@H](C)OC1=O'] 
    GG_i = Chem.MolFromSmiles(ring_smiles[0])
    LL_1 = Chem.MolFromSmiles(ring_smiles[1])
    LL_2 = Chem.MolFromSmiles(ring_smiles[2])
    LL_3 = Chem.MolFromSmiles(ring_smiles[3])   
    polymer = Chem.MolFromSmiles('C[C@H](O)C(=O)O[C@H](C)C(=O)O[I+]') if sequence[0:2]=='LL' else Chem.MolFromSmiles('[I+]OC(=O)COC(=O)CO')
    LA_count=0 if sequence[0:2]=='GG' else 2
    GA_count=0 if sequence[0:2]=='LL' else 2
    for i in range(0, len(sequence)-1,2):
        if sequence[i+2:i+4] == 'LL':
            polymer = reaction.RunReactants((polymer, LL_1))[0][0]
            Chem.SanitizeMol(polymer)
            LA_count+=2
        
        elif sequence[i+2:i+4] == 'GG':
            polymer = reaction.RunReactants((polymer, GG_i))[0][0]
            Chem.SanitizeMol(polymer)
            GA_count+=2

        elif sequence[i+2:i+4] == 'GL' or sequence[i+2:i+4] == 'LG':
            raise AttributeError("Check sequence. Input format is LLGG (blocks of dimers) and sequence entered is "+ sequence)
            
    
    LA_ratio = round((LA_count/len(sequence))*100,2)
    GA_ratio = round((GA_count/len(sequence))*100,2)
    if terminal == 'hydroxyl':
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('[OH]'), Chem.MolFromSmiles('O'))[0]
        Chem.AddHs(polymer)
    elif terminal == 'carboxyl':
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('[OH]'), Chem.MolFromSmiles('OC(=O)[OH]'))[0]
    elif terminal == 'ester':
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('[OH]'), Chem.MolFromSmiles('OC'))[0]
    else:
        raise ValueError("terminals accepts one of 3 arguments - hydroxyl, carboxyl or ester")
    polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('O[I]'), Chem.MolFromSmarts('[O]'))[0]
    Chem.SanitizeMol(polymer)
    return polymer, LA_ratio, GA_ratio




def blockiness_PLGA(sequence):  
    """
    Calculate the blockiness and average block length of a PLGA sequence.

    This function calculates the blockiness of a PLGA sequence by counting the occurrences of 'GG' and 'GL' or 'LG' in the sequence.
    It also calculates the average block length of 'G' and 'L' in the sequence.

    Parameters
    ----------
    sequence : str
        A string representing the PLGA sequence. 'G' represents one type of monomer and 'L' represents another type.

    Returns
    -------
    tuple
        A tuple containing the following elements:
        - blockiness (float): The blockiness of the PLGA sequence. Calculated as the ratio of 'GG' to 'GL' or 'LG'.
        - block_length_G (float): The average block length of 'G' in the sequence.
        - block_length_L (float): The average block length of 'L' in the sequence.

    Notes
    -----
    If the sequence does not contain both 'G' and 'L', the function returns a string indicating that the molecule is not a co-polymer.
    """

    if 'G' in sequence and 'L' in sequence:
        LG = sequence.count('LG')
        GG = sequence.count('GG')
        GL = sequence.count('GL')
        LL = sequence.count('LL')
        if 'GL' in sequence:
            blockiness = GG/GL
        else:
            blockiness = GG/LG
        
        block_list_G = [x for x in sequence.split('L') if x!='']
        block_length_G = mean([len(b) for b in block_list_G])
        
        block_list_L = [x for x in sequence.split('G') if x!='']
        block_length_L = mean([len(b) for b in block_list_L])
        return blockiness, block_length_G, block_length_L

    else:
        return 'Molecule is not a co-polymer, no blockiness calculation performed', 0, len(sequence)

def calculate_box_components_PLGA(chains, sequence, monomers = ['OC(=O)CO', 'C[C@@H](C(=O)[OH])O'], salt_concentration = 0.1 * unit.mole / unit.liter, residual_monomer = 0.00):
    """
    Calculate the components required to construct a simulation box for a given set of molecular chains.

    This function considers the salt concentration and residual monomer concentration to determine the quantity of each molecule type required.
    It is adapted from the OpenFF Toolkit Packmol wrapper's solvate_topology function.

    Parameters
    ----------
    chains : list
        A list of molecular chains to be included in the simulation box.
    sequence : str
        A string representing the sequence of the molecular chains. 'G' and 'L' represent different types of monomers.
    salt_concentration : float, optional
        The desired salt concentration in the simulation box. Defaults to 0.1 M.
    residual_monomer : float, optional
        The desired residual monomer concentration in the simulation box. Defaults to 0.05.

    Returns
    -------
    tuple
        A tuple containing the following elements:
        - molecules (list): A list of molecules to be included in the simulation box.
        - number_of_copies (list): A list indicating the quantity of each molecule to be included in the simulation box.
        - topology (openff.toolkit.topology.Topology): The topology of the simulation box.
        - box_vectors (numpy.ndarray): The vectors defining the dimensions of the simulation box.

    Notes
    -----
    This function is adapted from the OpenFF Toolkit Packmol wrapper's solvate_topology function.
    """
    # Function implementation here
    from openff.toolkit.topology import Molecule, Topology
    from openff.interchange.components._packmol import UNIT_CUBE, pack_box, RHOMBIC_DODECAHEDRON, solvate_topology
    from openff.interchange.components._packmol import _max_dist_between_points, _compute_brick_from_box_vectors, _center_topology_at
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
    tolerance= 2.0 * unit.nanometer
    
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
    water_to_add = int(round((solvent_mass - nacl_mass) / water_mass).m_as(unit.dimensionless).round())
    
    # Neutralise the system by adding and removing salt
    solute_charge = sum([molecule.total_charge for molecule in topology.molecules])
    na_to_add = int(round(np.ceil(nacl_to_add - solute_charge.m / 2.0)))
    cl_to_add = int(round(np.floor(nacl_to_add + solute_charge.m / 2.0)))
    
    rolling_mass=0
    for m in topology.molecules:
        rolling_mass += sum(atom.mass for atom in m.atoms)
    rolling_mass += nacl_mass * nacl_to_add
    rolling_mass += water_mass * water_to_add
    
    # residual monomer to add
    mass_to_add = (rolling_mass.magnitude/1-residual_monomer) * residual_monomer
    
    
    if 'L' in sequence and 'G' in sequence:
        A = Molecule.from_smiles(monomers[0])
        A_mass = sum([atom.mass for atom in A.atoms])
        B = Molecule.from_smiles(monomers[1])
        B_mass = sum([atom.mass for atom in B.atoms])
        for r in range(0,100):
            if (r * A_mass.magnitude) + (r * B_mass.magnitude) <= mass_to_add:
                A_to_add = r
                B_to_add = r
            else:
                break
        residual_monomer_actual = ((A_to_add * A_mass.magnitude + B_to_add * B_mass.magnitude) / rolling_mass.magnitude)
        molecules = [water, na, cl, A, B]
        number_of_copies=[water_to_add, na_to_add, cl_to_add, A_to_add, B_to_add]
    
    elif 'L' in sequence and 'G' not in sequence:
        A = Molecule.from_smiles(monomers[0])
        A_mass = sum([atom.mass for atom in A.atoms])
        B = Molecule.from_smiles('C')
        for r in range(0,100):
            if r * A_mass.magnitude <= mass_to_add:
                A_to_add = r
            else:
                break
        B_to_add = 0
        residual_monomer_actual = ((A_to_add * A_mass.magnitude) / rolling_mass.magnitude)
        molecules = [water, na, cl, A, B]
        number_of_copies=[water_to_add, na_to_add, cl_to_add, A_to_add, B_to_add]
    

    return molecules, number_of_copies, topology, box_vectors, residual_monomer_actual



#Class object for PLGA system
class PLGA_system:
    from openeye import oechem
    from openff.units import unit
    from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
    from functools import reduce
    from statistics import mean
    from rdkit.Chem.Descriptors import ExactMolWt
    from openff.interchange import Interchange
    from openff.interchange.components._packmol import UNIT_CUBE, pack_box
    from swiftpol.demo import build_PLGA_ring, blockiness_PLGA, calculate_box_components_PLGA
    from swiftpol.build import PDI
    from swiftpol import build
    from rdkit.Chem import AllChem


    gen_rxn = AllChem.ReactionFromSmarts('[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]')
    def __init__(self, perc_lactide_target, length_target, blockiness_target, terminals, num_chains):
        """
        A class used to represent a poly-lactide-(co)-glycolide (PLGA) polymer chain system.

        Attributes
        ----------
        lactide_target : float
            The target percentage of lactide in the polymer.
        length_target : int
            The target length of the polymer chain.
        blockiness_target : float
            The target blockiness of the polymer chain.
        terminals : str
            The end groups of the polymer.
        sequence : str
            The sequence of the polymer chain.
        chains : list
            A list of molecular chains in the system.
        chain_rdkit : list
            A list of RDKit molecule objects representing the chains.
        mol_weight_average : float
            The average molecular weight of the chains.
        PDI : float
            The Polydispersity Index of the chains.
        Mn : float
            The number-average molecular weight of the chains.
        Mw : float
            The weight-average molecular weight of the chains.
        num_chains : int
            The number of chains in the system.
        perc_lactide_actual : list
            A list of the actual percentage of lactide in each chain.
        lactide_actual : float
            The actual average percentage of lactide in the chains.
        length_average : float
            The average length of the chains.
        lengths : list
            A list of the lengths of the chains.
        min_length : int
            The minimum length of the chains.
        max_length : int
            The maximum length of the chains.
        blockiness_list : list
            A list of the blockiness of each chain.
        mean_blockiness : float
            The average blockiness of the chains.
        G_block_length : float
            The average block length for 'G' in the chains.
        L_block_length : float
            The average block length for 'L' in the chains.
        """
        self.lactide_target = perc_lactide_target
        self.length_target = length_target
        self.blockiness_target = blockiness_target
        self.terminals = terminals
        def spec(sequence, blockiness): #Define limits of lactide percentage and blockiness from input
            actual_lac = (sequence.count('L')/len(sequence))*100
            return actual_lac > perc_lactide_target*0.90 and actual_lac < perc_lactide_target*1.10 and blockiness>blockiness_target*0.90 and blockiness<blockiness_target*1.10
        chains = []
        perc_lactide_actual = []
        out_of_spec = 0
        chains_rdkit = []
        lengths = []
        blockiness_list = []
        GBL = []
        LBL = []
        self.monomers = ['OC(=O)CO', 'C[C@@H](C(=O)[OH])O']
        #First round of building
        for x in range(num_chains):
            length_actual = np.random.normal(length_target, 0.5)
            sequence = reduce(lambda x, y: x + y, np.random.choice(['LL', 'GG'], size=(int(length_actual/2),), p=[perc_lactide_target/100,1-(perc_lactide_target/100)]))
            blockiness = blockiness_PLGA(sequence)[0]
            if spec(sequence, blockiness)==True:
                reaction = build_PLGA_ring(sequence=sequence, terminal=terminals)
                lengths.append(int(length_actual))
                chains_rdkit.append(reaction[0])
                chain = Molecule.from_rdkit(reaction[0])
                chains.append(chain)
                perc_lactide_actual.append(reaction[1])
                blockiness_list.append(blockiness)
                GBL.append(blockiness_PLGA(sequence)[1])
                LBL.append(blockiness_PLGA(sequence)[2])
            else:
                out_of_spec +=1
        #Second round of building
        while out_of_spec >0:
            length_actual = np.random.normal(length_target, 0.5)
            sequence = reduce(lambda x, y: x + y, np.random.choice(['LL', 'GG'], size=(int(length_actual/2),), p=[perc_lactide_target/100,1-(perc_lactide_target/100)]))
            blockiness = blockiness_PLGA(sequence)[0]
            if spec(sequence, blockiness)==True:
                reaction = build_PLGA_ring(sequence=sequence, terminal=terminals)
                lengths.append(int(length_actual))
                chains_rdkit.append(reaction[0])
                chain = Molecule.from_rdkit(reaction[0])
                chains.append(chain)
                perc_lactide_actual.append(reaction[1])
                blockiness_list.append(blockiness)
                GBL.append(blockiness_PLGA(sequence)[1])
                LBL.append(blockiness_PLGA(sequence)[2])
                out_of_spec-=1
        self.sequence = sequence
        self.chains = chains
        self.chain_rdkit = chains_rdkit
        self.mol_weight_average = round(mean([ExactMolWt(c) for c in chains_rdkit]),2)
        self.PDI, self.Mn, self.Mw = PDI(chains_rdkit)
        self.num_chains = len(chains)
        self.perc_lactide_actual = perc_lactide_actual
        self.lactide_actual = mean(perc_lactide_actual)
        self.length_average = mean(lengths)
        self.lengths = lengths
        self.min_length = min(lengths)
        self.max_length = max(lengths)
        self.blockiness_list = blockiness_list
        self.mean_blockiness = mean(blockiness_list)
        self.G_block_length = mean(GBL)
        self.L_block_length = mean(LBL)
        print('System built!, size =', self.num_chains)

    def generate_conformers(self):
        """
        Generate conformers for each polymer chain in the system.

        This method uses the OpenFF toolkit to generate conformers for each polymer chain in the system. 
        It first checks if the OpenEye toolkit is licensed and available. If it is, it uses the OpenEyeToolkitWrapper 
        to generate conformers. Otherwise, it falls back to using the RDKitToolkitWrapper. Each chain is processed 
        to generate a single conformer, and unique atom names are assigned to each chain.

        Raises:
        ImportError: If neither RDKit nor OpenEye toolkits are available.
        """
        from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
        # Generate conformers using OpenFF toolkit wrapper
        for chain in self.chains:
            num = self.chains.index(chain)
            if oechem.OEChemIsLicensed():
                object = OpenEyeToolkitWrapper()
            else:
                object = RDKitToolkitWrapper()
            object.generate_conformers(molecule=chain, n_conformers=1)
            chain.generate_unique_atom_names()
            self.chains[num] = chain
    
    def charge_system(self):
        """
        Assign partial charges to each polymer chain in the system.

        This method uses the OpenFF NAGL toolkit to assign partial charges to each polymer chain in the system.
        It initializes a NAGLToolkitWrapper and assigns partial charges using the "openff-gnn-am1bcc-0.1.0-rc.2.pt" model.

        The method iterates over each chain in the `self.chains` list and assigns partial charges to the chain.

        Raises:
        ImportError: If the NAGL toolkit is not available.
        """
        from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        ntkw = NAGLToolkitWrapper()
        for chain in self.chains:
            ntkw.assign_partial_charges(chain, "openff-gnn-am1bcc-0.1.0-rc.2.pt")


    def solvate_system(self, resid_monomer, salt_concentration):
        """
        Build a solvated system using packmol functions.

        This method calculates the components needed to solvate the polymer system, including the molecules, 
        number of copies, topology, and box vectors. It then uses the `pack_box` function from the OpenFF 
        Interchange toolkit to create a solvated system.

        Parameters:
        resid_monomer (str): The residual monomer to be used in the system.
        salt_concentration (float): The concentration of salt to be added to the system.

        Returns:
        solvated_system: The solvated system generated by packmol.

        Raises:
        ImportError: If the OpenFF Interchange toolkit is not available.
        """
        from openff.interchange.components._packmol import pack_box
        
        molecules, number_of_copies, topology, box_vectors, resid_monomer_actual = calculate_box_components_PLGA(
            chains=self.chains,
            sequence=self.sequence, 
            residual_monomer=resid_monomer,
            salt_concentration=salt_concentration,
            monomers=self.monomers
        )
        self.residual_monomer = resid_monomer_actual
        self.solvent_comp = molecules
        self.num_copies_solvent = number_of_copies
        self.box_vectors = box_vectors
        solvated_system = pack_box(
            molecules=molecules,
            number_of_copies=number_of_copies,
            solute=topology,
            box_vectors=box_vectors,
            center_solute='BRICK'
        )
        return solvated_system
    

    def build_bulk(self, resid_monomer, salt_concentration=0 * unit.mole / unit.liter):
        """
        Build a bulk system using packmol functions.

        This method constructs a bulk system by packing the polymer chains into a box using the packmol algorithm.
        It calculates the topology from the polymer chains, determines the maximum distance between points in the 
        solute to set the box size, and then uses the `pack_box` function to create the bulk system.

        Parameters:
        resid_monomer (str): The residual monomer to be used in the system.
        salt_concentration (Quantity, optional): The concentration of salt to be added to the system. 
                                                   Defaults to 0 mole/liter.

        Returns:
        bulk_system: The bulk system generated by packmol.

        Raises:
        ImportError: If the OpenFF Interchange toolkit is not available.
        """

        solute_length = max(_max_dist_between_points(self.chains[i].to_topology().get_positions()) for i in range(len(self.chains)))
        box_vectors = UNIT_CUBE * solute_length
        bulk_system = pack_box(molecules=self.chains,
                                number_of_copies=[3 for i in range(len(self.chains))],
                                box_shape=UNIT_CUBE,
                                box_vectors=box_vectors,
                                center_solute='BRICK')
        
        return bulk_system
    
def build_starting_polymer(sequence='ABBBBBBBBBA', 
                           monomers=["IC=CC(=O)OI", "IOCCOI"], 
                           polymer_reaction=AllChem.ReactionFromSmarts("[C:1]-[O:2]-[I:3].[C:4]-[O:5]-[I:6]>>[C:1]-[O:2]-[C:4].[I:3]-[I:6].[O:5]"), 
                           activation_reaction=AllChem.ReactionFromSmarts('[C:1]=[C:2].[Cl:3].[I:4]>>[C:1](-[I:4])-[C:2]-[Cl:3]')):
    """
    Builds a polymer and adds reactive sites.

    Parameters:
        - sequence (str): Polymer sequence
        - monomers (list): List of monomers (SMILES format)
        - polymer_reaction (rdkit reaction): Reaction to form the polymer
        - activation_reaction (rdkit reaction): Reaction to introduce reactive sites

    Returns:
        - starting_polymer (rdkit Mol): PEGDA polymer with reactive sites
    """
    from swiftpol import build
    #Build PEGDA polymer
    polymer = build.build_polymer(sequence=sequence, monomer_list=monomers, reaction=polymer_reaction)

    activatorCl = Chem.MolFromSmiles('Cl')
    activatorI = Chem.MolFromSmiles('I')

    # Apply activation reaction twice
    starting_polymer = activation_reaction.RunReactants([polymer, activatorCl, activatorI])[0][0]
    Chem.SanitizeMol(starting_polymer)

    starting_polymer = activation_reaction.RunReactants([starting_polymer, activatorCl, activatorI])[0][0]
    Chem.SanitizeMol(starting_polymer)

    Chem.AddHs(starting_polymer)
    return starting_polymer

def replace_halogens_with_hydrogens(mol):
    """
    Replaces all halogen atoms (F, Cl, Br, I) in an RDKit molecule with hydrogen atoms.

    This function iterates through the atoms in the given molecule, identifies halogen atoms 
    (fluorine, chlorine, bromine, iodine), removes them, and replaces them with explicit hydrogens. 
    The resulting molecule is sanitized and updated to ensure chemical validity.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        The input RDKit molecule from which halogens will be replaced with hydrogens.

    Returns
    -------
    rdkit.Chem.Mol
        A new RDKit molecule with all halogens replaced by explicit hydrogens.

    Notes
    -----
    - The function uses an editable RDKit molecule (RWMol) to modify the structure.
    - The molecule is sanitized after halogen removal to ensure chemical correctness.
    - Explicit hydrogens are added to the resulting molecule.
    """
    editable = Chem.RWMol(mol)   #Creates an editable version of an RDKit molecule
    halogen_atomic_nums = [9, 17, 35, 53]  #F, Cl, Br, I

    for atom in list(editable.GetAtoms()):  #list() to avoid issues
        if atom.GetAtomicNum() in halogen_atomic_nums:
            idx = atom.GetIdx()
            #Remove the halogen
            editable.RemoveAtom(idx)

    #Sanitize and add hydrogens back
    mol_no_halogens = editable.GetMol()
    mol_no_halogens.UpdatePropertyCache()
    Chem.SanitizeMol(mol_no_halogens)

    #Add explicit hydrogens
    mol_with_h = Chem.AddHs(mol_no_halogens)
    return mol_with_h

def build_crosslinked_polymer(starting_polymer, 
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

    This is an experimental capability - proceed with caution.

    This function simulates the growth of a polymer network by either adding linear chains 
    or crosslinking existing polymer chains. The user can specify the number of iterations 
    and the probabilities of crosslinking versus linear chain addition. The resulting polymer 
    is returned with halogens replaced by hydrogens for chemical completeness.

    Parameters
    ----------
    starting_polymer : rdkit.Chem.Mol
        The initial polymer molecule to start the crosslinking or chain addition process.
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
    from swiftpol.demo import replace_halogens_with_hydrogens
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