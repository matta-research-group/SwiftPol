import pytest
from Parameterise import charge_polymer

def test_charge_polymer():
    # Test case 1: charge_scheme = 'AM1_BCC'
    polymer = 'CCO'
    charge_scheme = 'AM1_BCC'
    expected_charges = [0.0, 0.0, 0.0]
    assert charge_polymer(polymer, charge_scheme) == expected_charges

    # Test case 2: charge_scheme = 'espaloma'
    polymer = 'CCO'
    charge_scheme = 'espaloma'
    expected_charges = [0.0, 0.0, 0.0]
    assert charge_polymer(polymer, charge_scheme) == expected_charges

    # Test case 3: invalid charge_scheme
    polymer = 'CCO'
    charge_scheme = 'invalid'
    with pytest.raises(AttributeError):
        charge_polymer(polymer, charge_scheme)