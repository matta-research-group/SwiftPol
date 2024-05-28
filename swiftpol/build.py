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

#import mdtraj as md
import nglview
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



#Ring opening polymerisation
def build_PLGA_ring(sequence, 
                    reaction = AllChem.ReactionFromSmarts('[I:1][O:2].[I:3][C:4]>>[C:4][O:2].[I:3][I:1]'),
                    terminal='hydroxyl'):
    '''Build a PLGA co-polymer of specified sequence and return the sanitized polymer, specify monomer joining scheme using reaction SMARTS
    takes a list of up to 2 monomers to create a co-polymer.
    This function takes the cyclic esters lactide and glycolide as constituent monomers
    Inputs:
    reaction = Reaction SMARTS rdkit chemical reaction object specifying the joining of 2 iodinated compounds into an ester
    sequence = string with sequence (L for Lactide, G for glycolic acid). For this function, sequence must be assembled as blocks of 2 monomers
        e.g. LLGGLLGG
    monomer input is a RDkit.Chem Mol.object, not a SMILES string and must contain I-I bridging the point of esterification
    Outputs:
    PLGA macromolecule as RDkit.Chem.Mol
    Lactide ratio %
    Glycolide ratio %
    '''
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


#Class object for PLGA system - will be extended for other co-polymers and homopolymers


class PLGA_system:
    from openeye import oechem
    from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
    from functools import reduce
    from statistics import mean
    from rdkit.Chem.Descriptors import ExactMolWt
    from openff.interchange import Interchange
    from openff.interchange.components._packmol import UNIT_CUBE, pack_box
    from swiftpol.build import build_PLGA_ring
    """
    A class used to represent a poly-lactide-(co)-glycolide polymer chain system.

    ...

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

    Methods
    -------
    charge_system():
        Assigns partial charges to the chains and generates conformers.
    build_system(resid_monomer):
        Builds the system using packmol functions, containing a set amount of residual monomer.
    """

    gen_rxn = AllChem.ReactionFromSmarts('[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]')
    def __init__(self, perc_lactide_target, length_target, blockiness_target, terminals, num_chains): #Terminals will specify the end groups of the polymer (WIP)
        self.lactide_target = perc_lactide_target
        self.length_target = length_target
        self.blockiness_target = blockiness_target
        self.terminals = terminals
        def spec(sequence, blockiness): #Define limits of lactide percentage and blockiness from input
            actual_lac = (sequence.count('L')/len(sequence))*100
            return actual_lac > perc_lactide_target*0.95 and actual_lac < perc_lactide_target*1.05 and blockiness>blockiness_target*0.95 and blockiness<blockiness_target*1.05
        chains = []
        perc_lactide_actual = []
        out_of_spec = 0
        chains_rdkit = []
        lengths = []
        blockiness_list = []
        GBL = []
        LBL = []
        #First round of building
        for x in range(num_chains):
            length_actual = np.random.normal(length_target, 0.5)
            sequence = reduce(lambda x, y: x + y, np.random.choice(['LL', 'GG'], size=(int(length_actual/2),), p=[perc_lactide_target/100,1-(perc_lactide_target/100)]))
            blockiness = blockiness_calc(sequence)[0]
            if spec(sequence, blockiness)==True:
                reaction = build_PLGA_ring(sequence=sequence, terminal=terminals)
                lengths.append(int(length_actual))
                chains_rdkit.append(reaction[0])
                chain = Molecule.from_rdkit(reaction[0])
                chains.append(chain)
                perc_lactide_actual.append(reaction[1])
                blockiness_list.append(blockiness)
                GBL.append(blockiness_calc(sequence)[1])
                LBL.append(blockiness_calc(sequence)[2])
            else:
                out_of_spec +=1
        #Second round of building
        while out_of_spec >0:
            length_actual = np.random.normal(length_target, 0.5)
            sequence = reduce(lambda x, y: x + y, np.random.choice(['LL', 'GG'], size=(int(length_actual/2),), p=[perc_lactide_target/100,1-(perc_lactide_target/100)]))
            blockiness = blockiness_calc(sequence)[0]
            if spec(sequence, blockiness)==True:
                reaction = build_PLGA_ring(sequence=sequence, terminal=terminals)
                lengths.append(int(length_actual))
                chains_rdkit.append(reaction[0])
                chain = Molecule.from_rdkit(reaction[0])
                chains.append(chain)
                perc_lactide_actual.append(reaction[1])
                blockiness_list.append(blockiness)
                GBL.append(blockiness_calc(sequence)[1])
                LBL.append(blockiness_calc(sequence)[2])
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

    
    def charge_system(self): #Only supports GNN charging
        from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
        toolkit_registry = EspalomaChargeToolkitWrapper()
        for chain in self.chains:
            num = self.chains.index(chain)
            chain.assign_partial_charges('espaloma-am1bcc', toolkit_registry=toolkit_registry)
            #Generate conformers using OpenFF toolkit wrapper
            #if oechem.OEChemIsLicensed():
            #object = OpenEyeToolkitWrapper()
            #else:
            object = RDKitToolkitWrapper()
            chain.generate_unique_atom_names()


    def build_system(self, resid_monomer):
        '''Builds system using packmol functions'''
        self.residual_monomer = resid_monomer
        #IN DEVELOPMENT




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
    


def blockiness_calc(sequence):  
    """
    Calculates the blockiness and average block length of a co-polymer sequence.

    This function takes a sequence of co-polymers represented as 'G' and 'L'. It calculates the blockiness of the sequence, which is the ratio of 'GG' to 'GL' or 'LG', and the average block length for 'G' and 'L'. If the sequence does not contain both 'G' and 'L', it is not considered a co-polymer, and the blockiness is set to 1.0.

    Parameters:
    sequence (str): A string representing the sequence of co-polymers. 'G' and 'L' represent different types of monomers.

    Returns:
    tuple: A tuple containing the blockiness (float), the average block length for 'G' (float), and the average block length for 'L' (float).
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
        blockiness = 1.0
        block_length_G = len(sequence) if 'G' in sequence else 0
        block_length_L = len(sequence) if 'L' in sequence else 0
        print('Molecule is not a co-polymer, no blockiness calculation performed')
        return blockiness, block_length_G, block_length_L
