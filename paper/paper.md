--- 
title: 'SwiftPol: A Python package for building and parameterizing *in silico* polymer systems'
tags: 
  - Python 
  - polymer 
  - force field 
  - molecular dynamics 
  - polydispersity 
authors: 
  - name: Hannah N. Turney 
    orcid: 0009-0002-3298-0309 
    affiliation: 1 
  - name: Micaela Matta 
    orcid: 0000-0002-9852-3154
    affiliation: 1 
affiliations: 
 - name: Department of Chemistry, King’s College London Strand Campus, East Wing, 33-41 Surrey St, London, WC2R 2ND, United Kingdom 
   index: 1 
   ror: 0220mzb33

date:  DD Month YYYY 
bibliography: paper.bib 
 
 
--- 
 
# Summary 

A polymer sample contains a natural degree of variation in its structure and non-uniformity between its chains, which influences the bulk material properties of the sample. This innate heterogeneity is often disregarded in the *in silico* study of a polymer system, resulting in divergence from experiments. This paper presents SwiftPol, a user-guided Python software for the automated generation of polydisperse polymer ensembles which reproduce the heterogeneity observed in real materials. 

# Statement of need 
MD simulations of polymers are often performed with uniform idealized systems that do not capture the heterogeneity of their experimental counterparts. The result of this misalignment is non-convergence between MD-derived polymer properties and experimental data, and these MD simulations can overlook key components of polymer physics such as polydispersity and semi-crystallinity  [@schmid_understanding_2023; @triandafilidi_monodmolecular_2016]. Studies have demonstrated that bulk polymer material properties such as glass transition temperature, hydrophobicity, and inherent viscosity are highly sensitive to variations in polydispersity, making it essential to account for this heterogeneity to capture the true physics of polymer systems [@li_influence_2016; @wan_effect_2021; @ochi_influence_2021]. 
Polymer MD studies showcase an assortment of approaches to manually incorporate polydispersity into their polymer chain builds [@andrews_structure_2020; @kawagoe_construction_2019; @stipa_molecular_2021]. Although effective for their associated applications, these manual approaches are not universally applicable to different polymer chemistries or are performed using proprietary software.

Open-source software packages designed to build *in silico* polymer chains are focused on the design of polymers at the monomer and single-chain scale [@davel_parameterization_2024; @klein_hierarchical_2016; @santana-bonilla_modular_2023]. 

Packages that have the capability to build multi-chain systems such as Polyply [@grunewald_polyply_2022], RadonPy [@hayashi_radonpy_2022], and Polymer Structure Predictor [@sahu_polymer_2022], do not contain in-built functions to incorporate molecular weight diversity, and put the onus on the user to generate an ensemble of polydisperse chain sequences prior to using their chain building functions. Polymatic [@abbott_polymatic_2013], an algorithm to simulate polymerization, can theoretically create polydisperse systems through its random bond formation scheme, but does not have to ability to simultaneously control other key attributes of polymers such as blockiness and monomer ratio. The python tool Hoobas [@girard_hoobas_2019] builds polydisperse systems of randomized polymer chains for coarse-grained models. However, there is not currently a software package available that contains in-built functions to integrate multiple smaller-scale characteristics (monomer ratio, blockiness, degree of polymerization, chain terminal) into \textbf{atomistic} computational polymer models, whilst effectively capturing the heterogeneity and polydispersity of real-life samples. The development of SwiftPol was driven by the need to fill this gap in multi-scale building functionality of existing polymer building packages, to enable the simulation of realistic polymer models whilst allowing close user control of multiple system characteristics.


SwiftPol uses open-source Python libraries RDKit [@landrum_rdkitrdkit_2024], OpenFF-interchange [@thompson_openff_2024], and OpenFF-toolkit [@wagner_openforcefieldopenff-toolkit_2024] to promote reproducibility and portability.  We have ensured that SwiftPol objects can be seamlessly integrated into existing open-source software built for parameterization and simulation, to allow the user to select their preferred force field, topology format, and engine. RDKit, OpenFF-interchange and OpenFF-toolkit enable the export of SwiftPol polymer ensembles directly to simulation engines, and to a range of MD-compatible file formats, including .pdb, .top, .prmtop, and .json.

Here, we will detail the development of SwiftPol - a user-guided Python tool for building representative polymer ensembles, and subsequent studies to show its relevance and performance. 


# Package Overview 

The SwiftPol `build` module contains Python functions to build both single polymer chains and polydisperse polymer chain ensembles. 

SwiftPol takes as an input the simplified molecular-input line-entry system (SMILES)[@daylight_chemical_information_systems_inc_daylight_2022-1] string of all co-monomers, as well as values representing the target average properties of the ensemble: monomer % composition (for copolymers), length, number of chains, blockiness (for blocky copolymers), terminals, residual monomer. The user must define the reaction SMARTS [@daylight_chemical_information_systems_inc_daylight_2022] which describes the polymerization reaction associated with their polymer chemistry. 

As depicted in \autoref{Figure 1}, SwiftPol generates an initial polymer chain with a chain length drawn from a normal distribution centered around the specified target length, along with a terminal group that corresponds to the chosen input. In the case of a block copolymer, a probability function is used to determine the ratio of monomers in the chain and the chain is passed to a second function which tests whether the values for blockiness and % monomer are within 10% of the input variable by default. The +/- 10% acceptance margin introduces polydispersity into the ensemble by ensuring a certain level of non-uniformity between polymer chains, without straying too far from the input value. The acceptance margin can be adjusted by the user to control build stringency and polydispersity in the SwiftPol ensemble.

If all tests are passed, the chain is appended to the Python polymer ensemble build object, and the associated properties of the chain are calculated and added as ensemble attributes. Otherwise, the chain is discarded, and the process is repeated. Once the ensemble size is satisfied, average properties are calculated using built-in SwiftPol functions.

![Flowchart showing the process of building a polymer ensemble using SwiftPol. Created in BioRender. Matta, M. (2024) https://BioRender.com/o66z317.\label{Figure 1}](Fig_1_Swiftpol.png) 

This approach allows for the generation of a polydisperse chain ensemble, meaning each chain displays different properties but the ensemble matches the target properties and distribution, as is observed in experimental polymer samples. 

SwiftPol also contains functions to generate conformers using RDKit [@landrum_rdkitrdkit_2024] or OpenEye Omega(license-dependent, academic license provided by OpenEye, Cadence Molecular Sciences)[@hawkins_conformer_2010;@openeye_cadence_molecular_sciences_omega_2025], and assign force field parameters to the polydisperse ensembles using the openff-interchange infrastructure. The user can export the chain ensemble to existing tools such as packmol [@martinez_packmol_2009] and PolyPly [@grunewald_polyply_2022] to generate initial positions for molecular dynamics, the latter of which is particularly well-suited for building amorphous configurations for amorphous multi-component systems.

# Application: building a poly(lactide-co-glycolide) ensemble

Using SwiftPol, we have successfully constructed polydisperse ensembles of poly(lactide-co-glycolide) (PLGA), a widely used biodegradable polymer. We used the molecular structures and properties of experimental PLGA products as input for SwiftPol building functions to create representative PLGA systems to be used for molecular dynamics simulations. By integrating experimental data, such as chain terminals, copolymer ratios of lactic and glycolic acid, and blockiness, we have been able to replicate the bulk characteristics of various commercial polymer products, namely polydispersity. 
A full example implementation of SwiftPol for building PLGA systems can be found in the [building a PLGA system example notebook.](Example_Notebooks/PLGA_demo.ipynb)
We used SwiftPol to build ‘product X’, a commercially available 75:25 LA:GA ester-terminated PLGA. Following the chain build, another SwiftPol function was used to calculate the appropriate box size for the unit cell, number of water molecules, NaCl ions, and residual monomer molecules to include in the simulation of a complete condensed polymer ensemble.
The input values for the SwiftPol builder, seen in \autoref{tab:Table 1}, were taken from quality assurance documents provided by the manufacturer of product X, except the value for blockiness which was measured experimentally by Sun et al [@sun_characterization_2022]. The system attributes assigned by SwiftPol to the completed condensed PLGA unit cell are in seen in \autoref{tab:Table 2}.



\begin{flushleft}
\begin{table}[h!]
\captionsetup{justification=raggedright,singlelinecheck=false}
\caption[Input parameters for SwiftPol PLGA builder function, for the building of product X.\label{tab:Table 1}]
\begin{tabular}{|l|l|}
\hline
\textbf{INPUT} & \textbf{VALUE} \\
\hline
\hline
SYSTEM SIZE & 3 \\
TARGET LACTIDE PROPORTION (\%) & 75 \\
DEGREE OF POLYMERIZATION (MONOMER) & 50 \\
TARGET CHAIN BLOCKINESS & 1.7 \\
TERMINAL & Ester \\
RESIDUAL MONOMER (\% W/W) & 0.05 \\
NACL CONCENTRATION (M) & 0.1 \\
\hline
\end{tabular}
\end{table}
\end{flushleft}



\begin{flushleft}
\begin{table}[h!]
\captionsetup{justification=raggedright,singlelinecheck=false}
\caption[SwiftPol system build attributes. \( \bar{x}_n \) = mean value of attribute across n chains.\label{tab:Table 2}]
\begin{tabular}{|l|l|}
\hline
\textbf{ATTRIBUTE} & \textbf{\( \bar{x}_n \)} \\
\hline
\hline
SYSTEM SIZE (CHAINS) & 3 \\
ACTUAL LACTIDE PROPORTION (\%) & 68.9 \\
AVERAGE CHAIN BLOCKINESS & 1.65 \\
AVERAGE MOLECULE WEIGHT (DALTON) & 3370 \\
AVERAGE CHAIN LENGTH (MONOMERS) & 50 \\
POLYDISPERSITY INDEX  & 1.68 \\
BUILD TIME (S)  & 1.4 \\
\hline
\end{tabular}
\end{table}
\end{flushleft}

# Performance Benchmarking
We determined whether SwiftPol can build polymer ensembles and chains with sizes that are relevant to the system scales of interest by performing a stress test. \autoref{Figure 2} shows measurements of the performance benchmarking results, illustrating that SwiftPol can build large-scale systems in a realistic time frame.


![A) Time, t, taken to build systems with a single-chain, ranging from a 10-mer to a 1000-mer. B) Time, t, taken to 50-mer chain build systems ranging from 10 chains to 250 chains.\label{Figure 2}](Fig_2_Swiftpol.png) 


# Conclusion 

We presented SwiftPol, an open-source Python package for building polydisperse *in silico* polymer ensembles. SwiftPol recreates core characteristics of bulk polymer materials like polydispersity, enabling the simulation of representative systems that capture key components of polymer physics. We have shown that building longer chains and larger systems, exceeding what would be appropriate for atomistic MD simulations, will not create a time bottleneck in the MD workflow. SwiftPol is a robust and scalable tool for the guided generation of polydisperse polymer mixtures, which can be easily integrated into existing open-source MD software, such as the OpenFF toolkit.

In future releases, we will expand SwiftPol to include options to control tacticity, provide in-package integration into widely-established packing software and offer a broader selection of solvation buffers.

# Defining Polymer Properties 

SwiftPol uses the following expressions to define key polymer properties. 

`Monomer ratio`, *R~m~*, is the ratio of monomer A to monomer B in an AB copolymer, shown in \autoref{equation 1} 

\begin{equation}\label{equation 1} 
\mathit{R_{m}} = \frac{n(A)}{n(A+B)} 
\end{equation} 

`Degree of polymerization`, *DOP*, is the mean polymer chain length in the system, shown in \autoref{equation 2}. 

\begin{equation}\label{equation 2} 
DOP = \overline{x}(nA+nB) 
\end{equation} 

`Number of chains`, *n~chains~*, is the total number of chains built by SwiftPol and appended to the object, shown in \autoref{equation 3}. 

\begin{equation}\label{equation 3} 
n_{chains} = \mbox{total number of chains built}
\end{equation} 

`Blockiness`, *b*, is a measurement of the distribution of monomers in an AB copolymer, shown in \autoref{equation 4}. 

\begin{equation}\label{equation 4} 
\mathit{b} = \frac{nB-B\,\mbox{bonds}}{nA-B\,\mbox{bonds}} 
\end{equation} 

`Residual monomer`, *M~resid~*, is the % of residual monomer molecules in the system, shown in \autoref{equation 5}. 

\begin{equation}\label{equation 5} 
\mathit{M_{resid}} = \frac{M_{w}(\mbox{resid})}{M_{w}(\mbox{carbon-containing compounds})} 
\end{equation} 

# User guidance: Dependencies

[YAML environment file listing all required packages to use SwiftPol](https://github.com/matta-research-group/SwiftPol/blob/main/Dev_tools/swiftpol.yml)

To use SwiftPol please download the following packages:

- RDKit

- openff-interchange

- openff-toolkit

- openff-nagl

- openeye-toolkits

- openff-units==0.2.2

- dgl==2.0.0


optional dependencies:

- espaloma-charge


for quick dependency and swiftpol install using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), run the following commands in bash
```bash
conda create -n swiftpol python=3.10 rdkit openff-interchange openff-toolkit openff-nagl dgl=2.0.0 openeye-toolkits openff-units=0.2.2 nglview -c conda-forge -c dglteam -c openeye
conda activate swiftpol
git clone https://github.com/matta-research-group/SwiftPol
cd SwiftPol
pip install -e .
```

# Acknowledgements 

Hannah Turney is supported by funding contributions from the United Kingdom Research and Innovation Biotechnology and Biological Sciences Research Council (grant ref. BB/T008709/1) and Johnson&Johnson Innovative Medicine. 

We acknowledge contributions and feedback from Jeffrey Wagner at the Open Force field consortium and Anusha Lalitha, David Hahn, and Gary Tresadern at Johnson&Johnson Innovative Medicine.

We acknowledge the use of King’s College London e-research Computational Research, Engineering and Technology Environment (CREATE) high-performance computing facility in the development and testing of SwiftPol [@kings_college_london_kings_2024].

We acknowledge the use of a free academic license provided by OpenEye Scientific, Santa Fe, NM [@hawkins_conformer_2010; @openeye_cadence_molecular_sciences_omega_2025].

# References 
