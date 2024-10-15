import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import espaloma_charge as espcharge
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
from openff import toolkit
from openff.toolkit.topology import Molecule
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper


def charge_polymer(polymer, charge_scheme):
    '''
    Calculate and return the partial charges of a polymer chain based on the specified charge scheme.

    Parameters:
    polymer: A polymer chain for which the charges are to be calculated.
    charge_scheme: A string that specifies the charge scheme to be used. It can be either 'AM1_BCC', 'espaloma', or 'NAGL'.

    Returns:
    The partial charges of the polymer chain according to the specified charge scheme.

    Raises:
    AttributeError: If the charge_scheme input is not 'AM1_BCC', 'NAGL', or 'espaloma'.
    '''
    
    if charge_scheme == 'AM1_BCC':
        openff_chain = Molecule.from_rdkit(polymer)
        openff_chain.generate_conformers()
        openff_chain.assign_partial_charges("am1bcc")
        return openff_chain.partial_charges
    elif charge_scheme == 'espaloma':
        chain_h = Chem.AddHs(polymer)
        return espcharge.charge(chain_h)
    elif charge_scheme == 'NAGL' and toolkit.__version__ >= '0.16.0':
        chain_h = Chem.AddHs(polymer)
        openff_chain = Molecule.from_rdkit(chain_h)
        ntkw = NAGLToolkitWrapper()
        ntkw.assign_partial_charges(openff_chain, "openff-gnn-am1bcc-0.1.0-rc.2.pt")
        return openff_chain.partial_charges.magnitude
    
    else:
        raise AttributeError("This function takes either 'AM1_BCC', 'NAGL', or 'espaloma' as charge_scheme input")
    
    


def forcefield_with_charge_handler(molecule, charge_method, forcefield = "openff-2.2.0.offxml", ensemble=False):
    '''
    Create a forcefield with a charge handler for a given molecule and charge method. The function can handle both individual molecules and ensembles of molecules.

    Parameters:
    molecule: An RDKit molecule object or SwiftPol Ensemble for which the forcefield is to be created.
    charge_method: A string that specifies the charge method to be used for the molecule.
    forcefield: A string that specifies the OpenFF forcefield to be used. Default is "openff-2.2.0.offxml". IF a non-OpenFF or bespoke force field is being used,
    the user can specify the path to the force field file (format = SMIRNOFF XML).
    ensemble: A boolean that specifies whether the input molecule is an ensemble of molecules. Default is False.

    Returns:
    An OpenFF ForceField object with the specified molecule's charges added to the LibraryCharges parameter.

    Raises:
    Exception: If the charge method is not supported.
    '''

    
    import numpy as np
    from openmm import unit
    from openff.toolkit.topology import Molecule
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.toolkit.typing.engines.smirnoff.parameters import LibraryChargeHandler
    from swiftpol.parameterize import charge_polymer
    
    if ensemble==False:
        openff_molecule = Molecule.from_rdkit(molecule)
        charges = charge_polymer(molecule, charge_method)
        openff_molecule.partial_charges = charges * unit.elementary_charge

        #Create library charge type from openff_molecule
        library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(openff_molecule)

        # Pull base force field
        forcefield = ForceField(forcefield) #If forcefield is a string, it will be treated as a file path 
        # Add charges to OpenFF force field
        forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)
        
    elif ensemble==True:
        # Pull base force field
        forcefield = ForceField(forcefield)
        #Remove duplicate chains
        chains = set(molecule.chains)
        for i in chains:
            charges = charge_polymer(i.to_rdkit(), charge_method)
            i.partial_charges = charges * unit.elementary_charge
        
            # Add charges to OpenFF force field
            library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(i)

            forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)
    
    return forcefield

def forcefield_with_residualmonomer(system, charge_method, forcefield):
    """
    Add charges to a forcefield for a system with residual monomers.

    This function calculates the charges for lactide and glycolide monomers using the specified charge method, 
    and adds these charges to the provided forcefield. If the system does not contain any residual monomers, 
    the function raises an AttributeError.

    Parameters:
    system (object): An object representing the system. The system should have a 'residual_monomer' attribute 
                     indicating the number of residual monomers in the system.
    charge_method (str): The method to use for calculating charges. This should be a string recognized by the 
                         'charge_polymer' function.
    forcefield (ForceField): The forcefield to which to add the charges, OpenFF ForceField() object.

    Returns:
    ForceField: The forcefield with the added charges.

    Raises:
    AttributeError: If the system does not contain any residual monomers.
    """
    import numpy as np
    from openmm import unit
    from openff.toolkit.topology import Molecule
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.toolkit.typing.engines.smirnoff.parameters import LibraryChargeHandler
    from swiftpol.parameterize import charge_polymer
    if system.residual_monomer > 0:
        charges_lac = charge_polymer(Chem.MolFromSmiles('C[C@@H](C(=O)[OH])O'), charge_method)
        charges_gly = charge_polymer(Chem.MolFromSmiles('OC(=O)CO'), charge_method)
        lac = Molecule.from_rdkit(Chem.MolFromSmiles('C[C@@H](C(=O)[OH])O'))
        gly = Molecule.from_rdkit(Chem.MolFromSmiles('OC(=O)CO'))
        lac.partial_charges = charges_lac * unit.elementary_charge
        gly.partial_charges = charges_gly * unit.elementary_charge
            
        # Add charges to OpenFF force field
        library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(lac)
        forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)
        library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(gly)
        forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)
        
    else:
        raise AttributeError("The system does not contain any residual monomer")

    return forcefield