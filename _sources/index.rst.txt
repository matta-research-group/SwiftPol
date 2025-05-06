.. swiftpol documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SwiftPol's documentation!
=========================================================

About SwiftPol
--------------

A polymer sample contains a natural degree of variation in its structure and non-uniformity between its chains, which is often disregarded in MD systems.

SwiftPol uses a statistical approach to build polydisperse polymer systems, allowing for a more realistic representation of polymer systems in MD.

SwiftPol.build.polymer_system takes as an input: 

- The simplified molecular-input line-entry system (SMILES) string of all co-monomers.

- Values of the target average properties of the ensemble: monomer % composition (for copolymers), length, number of chains, blockiness (for blocky copolymers), terminals, residual monomer. 

- Reaction SMARTS which describes the polymerization reaction associated with their polymer chemistry.

...and build a polydisperse chain ensemble.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   api



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
