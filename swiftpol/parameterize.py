import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from openff import toolkit
from openff.toolkit.topology import Molecule


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
     from openff import toolkit
    from openff.toolkit.topology import Molecule
    if charge_scheme == 'AM1_BCC':
        openff_chain = Molecule.from_rdkit(polymer)
        openff_chain.generate_conformers()
        openff_chain.assign_partial_charges("am1bcc")
        return openff_chain.partial_charges
    elif charge_scheme == 'NAGL' and toolkit.__version__ < '0.16.0':
        raise ModuleNotFoundError("Installed version of openff-toolkit is below what is required to use NAGL. Please update to v.0.16.0")
    elif charge_scheme == 'NAGL' and toolkit.__version__ >= '0.16.0':
        try:
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        except:
            raise ImportError("The package openff-nagl is not installed. You will not be able to use NAGL.")
        chain_h = Chem.AddHs(polymer)
        openff_chain = Molecule.from_rdkit(chain_h)
        ntkw = NAGLToolkitWrapper()
        ntkw.assign_partial_charges(openff_chain, "openff-gnn-am1bcc-0.1.0-rc.2.pt")
        return openff_chain.partial_charges.magnitude
    elif charge_scheme == 'espaloma':
        try:
            import espaloma_charge as espcharge
        except ImportError:
            raise ImportError("The package espaloma-charge is not installed. You will not be able to use EspalomaCharge.")
        chain_h = Chem.AddHs(polymer)
        return espcharge.charge(chain_h)

    
    else:
        raise AttributeError("This function takes either 'AM1_BCC', 'NAGL', or 'espaloma' as charge_scheme input")
    


    


def forcefield_with_charge_handler(molecule, charge_method, forcefield = "openff-2.2.0.offxml", ensemble=False):
    '''
    Create a forcefield with a charge handler for a given molecule and charge method. The function can handle both individual molecules and ensembles of molecules.

    Parameters:
    molecule: An RDKit molecule object or a list of RDKit molecule objects for which the forcefield is to be created.
    charge_method: A string that specifies the charge method to be used for the molecule.
    forcefield: A string that specifies the forcefield to be used. Default is "openff-2.2.0.offxml".
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
        
        # Add charges to OpenFF force field
        library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(openff_molecule)

        # Pull base force field
        forcefield = ForceField(forcefield)
        forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)
        
    elif ensemble==True:
        # Pull base force field
        forcefield = ForceField(forcefield)
        for i in molecule.chain_rdkit:
            openff_molecule = Molecule.from_rdkit(i)
            charges = charge_polymer(i, charge_method)
            openff_molecule.partial_charges = charges * unit.elementary_charge
        
            # Add charges to OpenFF force field
            library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(openff_molecule)

            forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)
    
    return forcefield