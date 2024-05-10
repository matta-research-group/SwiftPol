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
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper, AmberToolsToolkitWrapper
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



def build_PLGA_ring(reaction, sequence, terminal='hydroxyl'):
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
    polymer = Chem.MolFromSmiles('C[C@H](O)C(=O)O[C@H](C)C(=O)O[I+]') if sequence[0:2]=='LL' else Chem.MolFromSmiles('[I+]OC(=O)COC(=O)CO')
    LA_count=0 if sequence[0:2]=='GG' else 2
    GA_count=0 if sequence[0:2]=='LL' else 2
    for i in range(0, len(sequence)-1,2):
        if sequence[i+2:i+4] == 'LL':
            polymer = i_rxn.RunReactants((polymer, LL_1))[0][0]
            Chem.SanitizeMol(polymer)
            LA_count+=2
        
        elif sequence[i+2:i+4] == 'GG':
            polymer = i_rxn.RunReactants((polymer, GG_i))[0][0]
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
    polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('O[I]'), Chem.MolFromSmiles('O'))[0]
    Chem.SanitizeMol(polymer)
    return polymer, LA_ratio, GA_ratio


def build_PLGA_linear(sequence):
    '''
    Build a polymer of specified sequence monomers and return the sanitized polymer.
    
    This function uses reaction SMARTS to specify the monomer joining scheme. It takes a list of up to 2 monomers to create a co-polymer. 
    The monomer input is a RDkit.Chem Mol object, not a SMILES string. The function also calculates and returns the ratio of Lactic Acid (LA) 
    and Glycolic Acid (GA) in the polymer.

    Parameters:
    sequence (list): A list of monomers to create a co-polymer. Each monomer is represented by a character ('L' for Lactic Acid and 'G' for Glycolic Acid).

    Returns:
    tuple: A tuple containing the sanitized polymer (RDkit.Chem Mol object), the ratio of Lactic Acid (LA) in the polymer, and the ratio of Glycolic Acid (GA) in the polymer.
    '''
    #Import monomers
    monomer_smiles = ['OC(=O)CO', 'C[C@@H](C(=O)[OH])O']
    glycolic = Chem.MolFromSmiles(monomer_smiles[0])
    lactate = Chem.MolFromSmiles(monomer_smiles[1])
    gen_rxn = AllChem.ReactionFromSmarts('[C:1][HO:2].[HO:3][C:4]>>[C:1][O:2][C:4].[O:3]')
    LA_count=0 if sequence[0]=='G' else 1
    GA_count=0 if sequence[0]=='L' else 1
    for i in range(len(sequence)-1):
        if sequence[i+1] == 'L':
            polymer = gen_rxn.RunReactants((polymer, lactate))[0][0]
            Chem.SanitizeMol(polymer)
            LA_count+=1
        
        elif sequence[i+1] == 'G':
            polymer = gen_rxn.RunReactants((polymer, glycolic))[0][0]
            Chem.SanitizeMol(polymer)
            GA_count+=1
    LA_ratio = round((LA_count/len(sequence))*100,2)
    GA_ratio = round((GA_count/len(sequence))*100,2)
    polymer = Chem.ReplaceSubstructs(polymer, Chem.MolFromSmarts('O[I]'), Chem.MolFromSmiles('O'))[0]
    Chem.SanitizeMol(polymer)
    
    return polymer, LA_ratio, GA_ratio


#Tests
import unittest
from rdkit import Chem
from build import build_PLGA_ring, build_PLGA_linear

class TestBuildPLGA(unittest.TestCase):
    def test_build_PLGA_ring(self):
        sequence = 'LLGG'
        polymer, LA_ratio, GA_ratio = build_PLGA_ring(sequence)
        self.assertIsInstance(polymer, Chem.rdchem.Mol)
        self.assertEqual(LA_ratio, 50.0)
        self.assertEqual(GA_ratio, 50.0)

    def test_build_PLGA_linear(self):
        sequence = 'LGLG'
        polymer, LA_ratio, GA_ratio = build_PLGA_linear(sequence)
        self.assertIsInstance(polymer, Chem.rdchem.Mol)
        self.assertEqual(LA_ratio, 50.0)
        self.assertEqual(GA_ratio, 50.0)

if __name__ == '__main__':
    unittest.main()
    
    
#Class object for polymer system
from functools import reduce
from statistics import mean
from rdkit.Chem.Descriptors import ExactMolWt
from openff.interchange import Interchange
from openff.interchange.components._packmol import UNIT_CUBE, pack_box

class PLGA_system:
    '''A poly-lactide-(co)-glycolide polymer chain system class'''
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
                reaction = build_PLGA_ring(gen_rxn, sequence, terminals)
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
                reaction = build_PLGA_ring(gen_rxn, sequence, terminals)
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

    
    def charge_system(self):
        toolkit_registry = EspalomaChargeToolkitWrapper()
        for chain in self.chains:
            num = self.chains.index(chain)
            chain.assign_partial_charges('espaloma-am1bcc', toolkit_registry=toolkit_registry)
            #Generate conformers using OpenFF toolkit wrapper
            object = OpenEyeToolkitWrapper()
            object.generate_conformers(molecule = chain, n_conformers=1)
            openff_chain.generate_unique_atom_names()


    def build_system(self, resid_monomer):
        '''Builds system using packmol functions'''
        self.residual_monomer = resid_monomer
        #IN DEVELOPMENT



def PDI(chains):
    #chains_rdkit = [Molecule.to_rdkit(chain) for chain in chains]
    mw_list = [ExactMolWt(chain) for chain in chains]  #_rdkit]
    #Mn = round(mean(mw_list),2)
    
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
        return 'Molecule is not a co-polymer, no blockiness calculation performed'
