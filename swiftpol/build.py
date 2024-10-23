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
from openeye import oechem

import openmm
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import get_data_file_path
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, AmberToolsToolkitWrapper
from openff.units import unit
from pandas import read_csv
import espaloma_charge as espcharge
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper

from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box


from functools import reduce
from statistics import mean
from rdkit.Chem.Descriptors import ExactMolWt
from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box

#Build polymer - generic
def build_polymer(sequence, monomer_list, reaction, terminal ='standard'):
    """
    Build a polymer using specified reaction sequence.

    This function takes a sequence of monomers, a list of corresponding monomer structures, and a reaction, 
    and builds a polymer using specified reaction sequence. The terminal group of the polymer can be 
    specified as 'hydroxyl' (default), 'carboxyl', or 'ester'.

    Parameters:
    sequence (str): A string representing the sequence of monomers. The sequence should be in the format 'AABBCCDD' 
                    (blocks of dimers).
    monomer_list (list): A list of SMILES strings representing the structures of the monomers. The order of the 
                         monomers in the list should correspond to the order of the monomers in the sequence.
    reaction (rdkit.Chem.rdChemReactions.ChemicalReaction): The reaction to use for the polymerization.
    terminal (str, optional): The terminal group of the polymer. Can be 'hydroxyl' (default), 'carboxyl', or 'ester'.

    Returns:
    rdkit.Chem.rdchem.Mol: The resulting polymer.

    Raises:
    AttributeError: If the sequence is not in the correct format.
    """
    monomers = {}
    for x in sorted(list(set(sequence))):
        ind = sorted(list(set(sequence))).index(x)
        monomers[x] = monomer_list[ind]
    polymer = Chem.MolFromSmiles('ICO')
    for i in range(0, len(sequence)):
        polymer = reaction.RunReactants((polymer, Chem.MolFromSmiles(monomers[sequence[i]])))[0][0]
        Chem.SanitizeMol(polymer)
            
    if terminal == 'hydroxyl':
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('[OH]'), Chem.MolFromSmiles('O'))[0]
        Chem.AddHs(polymer)
    elif terminal == 'carboxyl':
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('[OH]'), Chem.MolFromSmiles('OC(=O)[OH]'))[0]
    elif terminal == 'ester':
        polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('[OH]'), Chem.MolFromSmiles('OC'))[0]
    polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('OC[I]'), Chem.MolFromSmiles('O'))[0]
    Chem.SanitizeMol(polymer)
    return polymer


def build_linear_copolymer(sequence, 
                           monomer_a_smiles, 
                           monomer_b_smiles,
                           reaction = AllChem.ReactionFromSmarts('[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]')):
    """
    Constructs a linear co-polymer from the provided sequence of monomers.

    This function takes a sequence of monomers represented as 'A' and 'B', and the SMILES strings of two monomers. It constructs a co-polymer based on the sequence, using the provided reaction SMARTS for joining the monomers. The function returns the sanitized polymer and the percentage composition of each monomer in the polymer.

    Parameters:
    sequence (str): A string representing the sequence of monomers. 'A' represents monomer_a and 'B' represents monomer_b.
    monomer_a_smiles (str): The SMILES string of monomer A.
    monomer_b_smiles (str): The SMILES string of monomer B.
    reaction (rdkit.Chem.rdChemReactions.ChemicalReaction, optional): The reaction SMARTS used for joining the monomers. Defaults to '[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]', representing a condensation polymerisation.

    Returns:
    tuple: A tuple containing the sanitized polymer (rdkit.Chem.rdchem.Mol), the percentage composition of monomer A (float), and the percentage composition of monomer B (float).
    """

    polymer = Chem.MolFromSmiles('OC(=O)I')
    A_count=0
    B_count=0
    A = Chem.MolFromSmiles(monomer_a_smiles)
    B = Chem.MolFromSmiles(monomer_b_smiles)
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            polymer = reaction.RunReactants((polymer, A))[0][0]
            Chem.SanitizeMol(polymer)
            A_count+=1
        
        elif sequence[i] == 'B':
            polymer = reaction.RunReactants((polymer, B))[0][0]
            Chem.SanitizeMol(polymer)
            B_count+=1
    A_ratio = round((A_count/len(sequence))*100,2)
    B_ratio = round((B_count/len(sequence))*100,2)
    polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('OC(=O)I'), Chem.MolFromSmarts('[O]'))[0]
    Chem.SanitizeMol(polymer)
    return polymer, A_ratio, B_ratio



def PDI(chains):
    """
    Calculates the Polydispersity Index (PDI), number-average molecular weight (Mn), and weight-average molecular weight (Mw) of a list of chains.

    This function takes a list of molecular chains and calculates the PDI, which is the ratio of Mw to Mn. It also calculates Mn, which is the sum of the molecular weights of the chains divided by the number of chains, and Mw, which is the sum of the product of the weight fraction and molecular weight of each chain.

    Parameters:
    chains (list): A list of molecular chains. Each chain is represented as an RDkit molecule object.

    Returns:
    tuple: A tuple containing the PDI (float), Mn (float), and Mw (float).
    """

    mw_list = [ExactMolWt(chain) for chain in chains]  #_rdkit]
    list = [round(mass) for mass in mw_list]
    Mi = set(list)
    NiMi = []
    for i in Mi:
        Ni = list.count(i)
        NiMi.append(i*Ni)
    sigNiMi = sum(NiMi)
    Mn = sigNiMi/len(mw_list)
    wf = [z/sigNiMi for z in NiMi]
    WiMi = [wf[n]*NiMi[n] for n in range(len(wf))]
    Mw = sum(WiMi)
    PDI = Mw/Mn
    return PDI, Mn, Mw
    



def blockiness_gen(sequence):  
    """
    Calculate the blockiness and average block length of a co-polymer sequence.

    This function calculates the blockiness of a co-polymer sequence by counting the occurrences of 'GG' and 'GL' or 'LG' in the sequence. 
    It also calculates the average block length of 'G' and 'L' in the sequence.

    Parameters:
    sequence (str): A string representing the co-polymer sequence. 'G' represents one type of monomer and 'L' represents another type.

    Returns:
    blockiness (float): The blockiness of the co-polymer sequence. Calculated as the ratio of 'GG' to 'GL' or 'LG'.
    block_length_G (float): The average block length of 'G' in the sequence.
    block_length_L (float): The average block length of 'L' in the sequence.

    If the sequence does not contain both 'G' and 'L', the function returns a string indicating that the molecule is not a co-polymer.
    """

    if 'A' in sequence and 'B' in sequence:
        AB = sequence.count('AB')
        BB = sequence.count('BB')
        BA = sequence.count('BA')
        AA = sequence.count('AA')
        if 'BA' in sequence:
            blockiness = BB/BA
        else:
            blockiness = BB/AB
        
        block_list_B = [x for x in sequence.split('A') if x!='']
        block_length_B = mean([len(b) for b in block_list_B])
        
        block_list_A = [x for x in sequence.split('B') if x!='']
        block_length_A = mean([len(b) for b in block_list_A])
        return blockiness, block_length_B, block_length_A

    else:
        return 'Molecule is not a co-polymer, no blockiness calculation performed', 0, len(sequence)

def calculate_box_components(chains, sequence, salt_concentration = 0.1 * unit.mole / unit.liter, residual_monomer = 0.00):
    """
    ADAPTED FROM OPENFF TOOLKIT PACKMOL WRAPPER SOLVATE_TOPOLOGY FUNCTION
    Calculates the components required to construct a simulation box for a given set of molecular chains.
    Considers the salt concentration and residual monomer concentration to determine the quantity of each molecule type required.

    Parameters:
    chains (list): A list of molecular chains to be included in the simulation box.
    sequence (str): A string representing the sequence of the molecular chains. 'G' and 'L' represent different types of monomers.
    salt_concentration (float, optional): The desired salt concentration in the simulation box. Defaults to 0.1 M.
    residual_monomer (float, optional): The desired residual monomer concentration in the simulation box. Defaults to 0.05.

    Returns:
    tuple: A tuple containing the following elements:
        - molecules (list): A list of molecules to be included in the simulation box.
        - number_of_copies (list): A list indicating the quantity of each molecule to be included in the simulation box.
        - topology (openff.toolkit.topology.Topology): The topology of the simulation box.
        - box_vectors (numpy.ndarray): The vectors defining the dimensions of the simulation box.
    """
    from openff.toolkit.topology import Molecule, Topology
    from openff.interchange.components._packmol import UNIT_CUBE, pack_box, RHOMBIC_DODECAHEDRON, solvate_topology
    from openff.interchange.components._packmol import _max_dist_between_points, _compute_brick_from_box_vectors, _center_topology_at
    #Create molecules for the purpose of later topology assembly
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
    
    topology = Topology.from_molecules(chains)
    nacl_conc=salt_concentration
    padding= 0.1 * unit.nanometer
    box_shape= UNIT_CUBE
    target_density= 1.0 * unit.gram / unit.milliliter
    tolerance= 2.0 * unit.nanometer
    
    # Compute box vectors from the solute length and requested padding
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
    
    
    mass_to_add = (rolling_mass.magnitude/1-residual_monomer) * residual_monomer
    lac = Molecule.from_smiles('C[C@@H](C(=O)[OH])O')
    lac_mass = sum([atom.mass for atom in lac.atoms])
    gly = Molecule.from_smiles('OC(=O)CO')
    gly_mass = sum([atom.mass for atom in gly.atoms])
    
    
    if 'L' in sequence and 'G' in sequence:
        for r in range(0,100):
            if (r * lac_mass.magnitude) + (r * gly_mass.magnitude) <= mass_to_add:
                lac_to_add = r
                gly_to_add = r
            else:
                break
    
    elif 'L' in sequence and 'G' not in sequence:
        for r in range(0,100):
            if r * lac_mass.magnitude <= mass_to_add:
                lac_to_add = r
            else:
                break
        gly_to_add = 0
    
    residual_monomer_actual = ((lac_to_add * lac_mass.magnitude + gly_to_add * gly_mass.magnitude) / rolling_mass.magnitude)

    molecules = [water, na, cl, lac, gly]
    number_of_copies=[water_to_add, na_to_add, cl_to_add, lac_to_add, gly_to_add]
    return molecules, number_of_copies, topology, box_vectors, residual_monomer_actual



#Class object for generic polymer system

class polymer_system:
    from openeye import oechem
    from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
    from functools import reduce
    from statistics import mean
    from rdkit.Chem.Descriptors import ExactMolWt
    from openff.interchange import Interchange
    from openff.interchange.components._packmol import UNIT_CUBE, pack_box
    from swiftpol.build import build_polymer, PDI, blockiness_gen, calculate_box_components
    from openff.units import unit
    from rdkit.Chem import AllChem
    import numpy as np

    def __init__(self, monomer_list, reaction, length_target, num_chains, terminals='standard', perc_A_target=100, blockiness_target=1.0, copolymer=False):
        """
        Initialize the polymer system and build the polymer chains.

        Parameters:
        monomer_list (list): List of monomers to be used in the polymerization.
        reaction (str): The type of reaction to be used for polymerization.
        length_target (float): The target length of the polymer chains.
        num_chains (int): The number of polymer chains to be generated.
        terminals (str, optional): The type of terminal groups to be used. Default is 'standard'.
        perc_A_target (float, optional): The target percentage of monomer A in the copolymer. Default is 100.
        blockiness_target (float, optional): The target blockiness of the copolymer. Default is 1.0.
        copolymer (bool, optional): Flag to indicate if the system is a copolymer. Default is False.

        Attributes:
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
            self.blockiness_target = blockiness_target
            self.A_target = perc_A_target
            def spec(sequence, blockiness): #Define limits of A percentage and blockiness from input
                actual_A = (sequence.count('A')/len(sequence))*100
                blockiness = blockiness_gen(sequence)[0]
                return actual_A > perc_A_target*0.90 and actual_A < perc_A_target*1.10 and blockiness>blockiness_target*0.90 and blockiness<blockiness_target*1.10
            
            blockiness_list = []
            out_of_spec = 0
            BBL = []
            ABL = []
        chains = []
        chains_rdkit = []
        lengths = []
        

        #First round of building
        
        if copolymer==True:
            for n in range(num_chains):
                length_actual = np.random.normal(length_target, 0.5)
                sequence = reduce(lambda x, y: x + y, np.random.choice(['A', 'B'], size=(int(length_actual/2),), p=[perc_A_target/100,1-(perc_A_target/100)]))
                blockiness = blockiness_gen(sequence)[0]
                if spec(sequence, blockiness)==True:
                    pol = build_polymer(sequence=sequence, monomer_list = monomer_list, reaction = reaction, terminal=terminals)
                    lengths.append(int(length_actual))
                    chains_rdkit.append(pol)
                    chain = Molecule.from_rdkit(pol)
                    chains.append(chain)
                    perc_A_actual.append((sequence.count('A')/len(sequence))*100)
                    blockiness_list.append(blockiness)
                    BBL.append(blockiness_gen(sequence)[1])
                    ABL.append(blockiness_gen(sequence)[2])
                else:
                    out_of_spec +=1
                #Second round of building
                while out_of_spec >0:
                    length_actual = np.random.normal(length_target, 0.5)
                    sequence = reduce(lambda x, y: x + y, np.random.choice(['A', 'B'], size=(int(length_actual/2),), p=[perc_A_target/100,1-(perc_A_target/100)]))
                    blockiness = blockiness_gen(sequence)[0]
                    if spec(sequence, blockiness)==True:
                        pol = build_polymer(sequence=sequence, monomer_list = monomer_list, reaction = reaction, terminal=terminals)
                        lengths.append(int(length_actual))
                        chains_rdkit.append(pol)
                        chain = Molecule.from_rdkit(pol)
                        chains.append(chain)
                        perc_A_actual.append((sequence.count('A')/len(sequence))*100)
                        blockiness_list.append(blockiness)
                        BBL.append(blockiness_gen(sequence)[1])
                        ABL.append(blockiness_gen(sequence)[2])
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
                sequence = reduce(lambda x, y: x + y, np.random.choice(['A', 'B'], size=(int(length_actual/2),), p=[perc_A_target/100,1-(perc_A_target/100)]))
                pol = build_polymer(sequence=sequence, monomer_list = monomer_list, reaction = reaction, terminal=terminals)
                lengths.append(int(length_actual))
                chains_rdkit.append(pol)
                chain = Molecule.from_rdkit(pol)
                chains.append(chain)
                perc_A_actual.append((sequence.count('A')/len(sequence))*100)
        self.sequence = sequence
        self.chains = chains
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
        
        molecules, number_of_copies, topology, box_vectors, resid_monomer_actual = calculate_box_components(
            chains=self.chains,
            sequence=self.sequence, 
            residual_monomer=resid_monomer,
            salt_concentration=salt_concentration
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
        solute_length = max(_max_dist_between_points(sys.chains[i].to_topology().get_positions()) for i in range(len(sys.chains)))
        box_vectors = UNIT_CUBE * solute_length
        bulk_system = pack_box(molecules=sys.chains,
                                number_of_copies=[3 for i in range(len(sys.chains))],
                                box_shape=UNIT_CUBE,
                                box_vectors=box_vectors,
                                center_solute='BRICK')
        
        return bulk_system

