API Documentation
=================

.. autosummary::
   :toctree: autosummary

   swiftpol.canvas

swiftpol.build
==============

.. automodule:: swiftpol.build
   :members:
   :undoc-members:
   :show-inheritance:

Functions
---------

build_polymer
-------------

.. autofunction:: swiftpol.build.build_polymer

Example:
    Here is an example of how to use the `build_polymer` function:

.. code-block:: python
    from swiftpol import build
    from rdkit.Chem import AllChem

    # Build a short PEG chain
    polymer = build.build_polymer(
        'AAAA',
        ['IOCCOI'], 
        AllChem.ReactionFromSmarts('[C:1]-[O:2]-[I:3].[C:4]-[O:5]-[I:6]>>[C:1]-[O:2]-[C:4].[I:3][O:5][I:6]')
    )
    polymer


