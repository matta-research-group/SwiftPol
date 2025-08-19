import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from openff import toolkit
from openff.toolkit.topology import Molecule


def charge_polymer(polymer, charge_scheme):
    """
    Calculate and return the partial charges of a polymer chain based on the specified charge scheme.

    Parameters
    ----------
    polymer : rdkit.Chem.rdchem.Mol
        A polymer chain for which the charges are to be calculated.
    charge_scheme : str
        A string that specifies the charge scheme to be used. It can be either 'AM1_BCC', 'espaloma', or 'NAGL'.

    Returns
    -------
    numpy.ndarray
        The partial charges of the polymer chain according to the specified charge scheme.

    Raises
    ------
    AttributeError
        If the charge_scheme input is not 'AM1_BCC', 'NAGL', or 'espaloma'.
    """
    # Function implementation here
    from openff import toolkit
    from openff.toolkit.topology import Molecule

    if charge_scheme == "AM1_BCC":
        openff_chain = Molecule.from_rdkit(polymer)
        openff_chain.generate_conformers()
        openff_chain.assign_partial_charges("am1bcc")
        return openff_chain.partial_charges.magnitude
    elif charge_scheme == "NAGL" and toolkit.__version__ < "0.16.0":
        raise ModuleNotFoundError(
            "Installed version of openff-toolkit is below what is required to use NAGL. Please update to v.0.16.0"
        )
    elif charge_scheme == "NAGL" and toolkit.__version__ >= "0.16.0":
        try:
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        except:
            raise ImportError(
                "The package openff-nagl is not installed. You will not be able to use NAGL."
            )
        chain_h = Chem.AddHs(polymer)
        openff_chain = Molecule.from_rdkit(chain_h)
        ntkw = NAGLToolkitWrapper()
        ntkw.assign_partial_charges(openff_chain, "openff-gnn-am1bcc-0.1.0-rc.2.pt")
        return openff_chain.partial_charges.magnitude
    elif charge_scheme == "espaloma":
        try:
            from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
        except ImportError:
            raise ImportError(
                "The package espaloma-charge is not installed. You will not be able to use EspalomaCharge."
            )
        chain_h = Chem.AddHs(polymer)
        openff_chain = Molecule.from_rdkit(chain_h)
        etkw = EspalomaChargeToolkitWrapper()
        openff_chain.assign_partial_charges("espaloma-am1bcc", toolkit_registry=etkw)
        return openff_chain.partial_charges.magnitude

    else:
        raise AttributeError(
            "This function takes either 'AM1_BCC', 'NAGL', or 'espaloma' as charge_scheme input"
        )


def charge_openff_polymer(openff_chain, charge_scheme, overwrite=True):
    """
    Assign partial charges to a polymer chain based on the specified charge scheme.

    This function assigns partial charges to a polymer chain using one of the following charge schemes:
    'AM1_BCC', 'NAGL', or 'espaloma'. It can optionally overwrite existing charges.

    Parameters
    ----------
    openff_chain : openff.toolkit.topology.Molecule
        The original polymer chain to which partial charges will be assigned.
    charge_scheme : str
        The charge scheme to use for assigning partial charges. Options are 'AM1_BCC', 'NAGL', or 'espaloma'.
    overwrite : bool, optional
        Whether to overwrite existing partial charges. Default is True.

    Returns
    -------
    openff.toolkit.topology.molecule.PartialChargeArray
        The partial charges assigned to the polymer chain.

    Raises
    ------
    ModuleNotFoundError
        If the installed version of openff-toolkit is below what is required to use NAGL.
    ImportError
        If the required package for the specified charge scheme is not installed.
    AttributeError
        If the charge_scheme input is not 'AM1_BCC', 'NAGL', or 'espaloma'.
    """
    # Function implementation here
    from openff import toolkit
    from openff.toolkit.topology import Molecule

    chain_copy = openff_chain
    if charge_scheme == "AM1_BCC":
        if chain_copy.conformers is None:
            chain_copy.generate_conformers()
        if overwrite:
            openff_chain.assign_partial_charges("am1bcc")
            return openff_chain.partial_charges
        else:
            chain_copy.assign_partial_charges("am1bcc")
            return chain_copy.partial_charges
    elif charge_scheme == "NAGL" and toolkit.__version__ < "0.16.0":
        raise ModuleNotFoundError(
            "Installed version of openff-toolkit is below what is required to use NAGL. Please update to v.0.16.0"
        )
    elif charge_scheme == "NAGL" and toolkit.__version__ >= "0.16.0":
        try:
            from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        except:
            raise ImportError(
                "The package openff-nagl is not installed. You will not be able to use NAGL."
            )
        ntkw = NAGLToolkitWrapper()
        if overwrite:
            ntkw.assign_partial_charges(openff_chain, "openff-gnn-am1bcc-0.1.0-rc.2.pt")
            return openff_chain.partial_charges
        else:
            ntkw.assign_partial_charges(chain_copy, "openff-gnn-am1bcc-0.1.0-rc.2.pt")
            return chain_copy.partial_charges
    elif charge_scheme == "espaloma":
        try:
            from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
        except ImportError:
            raise ImportError(
                "The package espaloma-charge is not installed. You will not be able to use EspalomaCharge."
            )
        etkw = EspalomaChargeToolkitWrapper()
        if overwrite:
            openff_chain.assign_partial_charges(
                "espaloma-am1bcc", toolkit_registry=etkw
            )
            return openff_chain.partial_charges
        else:
            chain_copy.assign_partial_charges("espaloma-am1bcc", toolkit_registry=etkw)
            return chain_copy.partial_charges

    else:
        raise AttributeError(
            "This function takes either 'AM1_BCC', 'NAGL', or 'espaloma' as charge_scheme input"
        )


def forcefield_with_charge_handler(
    molecule, charge_method, forcefield="openff-2.2.0.offxml", ensemble=False
):
    """
    Create a forcefield with a charge handler for a given molecule and charge method.

    The function can handle both individual molecules and ensembles of molecules.

    Parameters
    ----------
    molecule : rdkit.Chem.rdchem.Mol or list of rdkit.Chem.rdchem.Mol
        An RDKit molecule object or a list of RDKit molecule objects for which the forcefield is to be created.
    charge_method : str
        A string that specifies the charge method to be used for the molecule.
    forcefield : str, optional
        A string that specifies the forcefield to be used. Default is "openff-2.2.0.offxml".
    ensemble : bool, optional
        A boolean that specifies whether the input molecule is an ensemble of molecules. Default is False.

    Returns
    -------
    openff.toolkit.typing.engines.smirnoff.ForceField
        An OpenFF ForceField object with the specified molecule's charges added to the LibraryCharges parameter.

    Raises
    ------
    Exception
        If the charge method is not supported.
    """
    # Function implementation here
    import numpy as np
    from openmm import unit
    from openff.toolkit.topology import Molecule
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.toolkit.typing.engines.smirnoff.parameters import LibraryChargeHandler
    from swiftpol.parameterize import charge_polymer

    if ensemble == False:
        openff_molecule = Molecule.from_rdkit(molecule)
        charges = charge_polymer(molecule, charge_method)
        openff_molecule.partial_charges = charges * unit.elementary_charge

        # Add charges to OpenFF force field
        library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(
            openff_molecule
        )

        # Pull base force field
        forcefield = ForceField(forcefield)
        forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)

    elif ensemble == True:
        # Pull base force field
        forcefield = ForceField(forcefield)
        for i in molecule.chain_rdkit:
            openff_molecule = Molecule.from_rdkit(i)
            charges = charge_polymer(i, charge_method)
            openff_molecule.partial_charges = charges * unit.elementary_charge

            # Add charges to OpenFF force field
            library_charge_type = LibraryChargeHandler.LibraryChargeType.from_molecule(
                openff_molecule
            )

            forcefield["LibraryCharges"].add_parameter(parameter=library_charge_type)

    return forcefield
