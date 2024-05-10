def charge_polymer(polymer, charge_scheme):
    '''return either AM1_BCC or espaloma charges for a polymer chain'''
    if charge_scheme == 'AM1_BCC':
        openff_chain = Molecule.from_rdkit(polymer)
        openff_polymer = openff_chain.generate_conformers()
        openff_polymer.assign_partial_charges("am1bcc")
        return openff_polymer.partial_charges
    elif charge_scheme == 'espaloma':
        chain_h = Chem.AddHs(polymer)
        return espcharge.charge(chain_h)
    else:
        raise AttributeError("This function takes either 'AM1_BCC' or 'espaloma' as charge_scheme input")