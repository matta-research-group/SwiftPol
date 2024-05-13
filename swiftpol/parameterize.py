import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import espaloma_charge as espcharge
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
from openff.toolkit.topology import Molecule

def print_polymer_charges(polymer, charge_scheme):
    '''
    Calculate and return the partial charges of a polymer chain based on the specified charge scheme.

    Parameters:
    polymer: A polymer chain for which the charges are to be calculated.
    charge_scheme: A string that specifies the charge scheme to be used. It can be either 'AM1_BCC' or 'espaloma'.

    Returns:
    The partial charges of the polymer chain according to the specified charge scheme.

    Raises:
    AttributeError: If the charge_scheme input is not 'AM1_BCC' or 'espaloma'.
    '''

    if charge_scheme == 'AM1_BCC':
        openff_chain = Molecule.from_rdkit(polymer)
        openff_chain.generate_conformers()
        openff_chain.assign_partial_charges("am1bcc")
        return openff_polymer.partial_charges
    elif charge_scheme == 'espaloma':
        chain_h = Chem.AddHs(polymer)
        return espcharge.charge(chain_h)
    else:
        raise AttributeError("This function takes either 'AM1_BCC' or 'espaloma' as charge_scheme input")
    

def assign_polymer_charges(polymer, charge_scheme):
    '''
    Assign either AM1_BCC or espaloma charges to a polymer chain.

    Parameters:
    polymer: A polymer chain for which the charges are to be assigned.
    charge_scheme: A string that specifies the charge scheme to be used. It can be either 'AM1_BCC' or 'espaloma'.

    Returns:
    None. The function modifies the polymer object in-place.

    Raises:
    AttributeError: If the charge_scheme input is not 'AM1_BCC' or 'espaloma'.
    '''
    
    if charge_scheme == 'AM1_BCC':
        openff_chain = Molecule.from_rdkit(polymer)
        openff_polymer = openff_chain.generate_conformers()
        openff_polymer.assign_partial_charges("am1bcc")
    elif charge_scheme == 'espaloma':
        chain_h = Chem.AddHs(polymer)
    else:
        raise AttributeError("This function takes either 'AM1_BCC' or 'espaloma' as charge_scheme input")