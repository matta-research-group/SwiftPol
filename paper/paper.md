--- 
title: ‘SwiftPol: A Python package for building and parameterizing *in silico* polymer systems' 
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
 - name: Department of Chemistry, King’s College London 
   index: 1 
   ror: ?? 
date: DD Month YYYY 
bibliography: paper.bib 
 
 
--- 
 
# Summary 
 

A polymer sample contains a natural degree of variation in its structure and non-uniformity between their chains, which influences the bulk material properties of the polymer. The innate heterogeneity is often disregarded in the *in silico* study of polymer systems. This paper presents ‘SwiftPol,’ a user-guided Python software for the automated generation of polydisperse polymer systems which are representative of experimental samples. 

We have ensured that SwiftPol can be seamlessly integrated into existing open-source software built for parameterization and simulation, to allow the user to use their preferred force field and engines. 

 
# Statement of need 
 
Polymer MD simulations are often performed with uniform idealized systems that do not capture the heterogeneity of their experimental counterparts. The result of this misalignment is non-convergence between MD-derived polymer properties and experimental data, and these MD simulations can overlook key components of polymer physics such as polydispersity and semi-crystallinity. Studies have demonstrated that polymer material properties are highly sensitive to variations in substructure, making it essential to account for this bulk diversity in order to capture the true physics of polymer systems `[@Wan:2001]`. 

Existing polymer MD studies showcase an assortment of approaches to manually incorporate polydispersity into their systems. However, there is not currently a software package available that effectively builds polymer systems that capture the naturally occurring heterogeneity and polydispersity. Here, we will detail our development of a user-guided python tool for building representative polymer systems, and subsequent studies to show its relevance and performance. 

SwiftPol uses open-source Python libraries ‘RDkit’ and ‘Open Force Field’ to promote accessibility.  

# Package Overview 

The SwiftPol build module contains Python functions to build both single polymer chains, and polydisperse polymer systems. 

SwiftPol takes as an input the simplified molecular-input line-entry system (SMILES) string of all monomer components in a polymer, as well as a panel of parameters representing the target average properties of the system; monomer % composition (for co-polymers), length, number of chains, blockiness (for blocky co-polymers), terminals, residual monomer. The user must define the reaction SMARTS which describes the polymerization reaction associated with their system. 

As depicted in Figure 1 \autoref{fig:1}, SwiftPol generates an initial polymer sequence that fits the chain length and terminal parameters exactly, using a probability function to determine the ratio of monomers in the chain. In the case of a block co-polymer, the chain is passed to a second function which tests whether the values for blockiness and % monomer are within 10% of the input variable. The +/- 10% acceptance margin introduces polydispersity into the system by ensuring a certain level of non-uniformity between polymer chains, without straying too far from the input value. 

If all tests are passed, the chain is appended to the Python polymer system object, and the associated actual properties of the chain are calculated and added as system attributes. Otherwise, the chain is discarded, and the process is repeated. Once the system size is satisfied, average system properties are calculated using in-built SwiftPol functions. 

![Figure 1. Flowchart showing the process of building a polymer system using SwiftPol.\label{fig:1}](Fig_1_swiftpol.png) 

This approach allows for the generation of a polydisperse system, meaning each chain displays the same average properties but the complete system exhibits the distributions observed in experimental polymer samples. 

SwiftPol also contains functions to generate conformers and assign force field parameters to the polydisperse systems. 

# Applications 

Using SwiftPol, we have successfully constructed polydisperse systems of poly(lactide-co-glycolide) (PLGA), a widely used biodegradable polymer. We used the molecular structures and properties of experimental PLGA products as input for SwiftPol building functions to create representative PLGA systems for molecular dynamics study. By integrating experimental data, such as chain terminals, copolymer ratios of lactic and glycolic acid, and blockiness, we have been able to replicate the bulk characteristics of various commercial polymer products, namely polydispersity. 

# Conclusion 

Considering key experimental properties of bulk polymer materials during *in silico* system building creates representative simulations which are more characteristic of experimental systems] 

 
# Defining Polymer Properties 

SwiftPol uses the following expressions to define key polymer properties. 

Monomer ratio, *R~m*, is the ratio of monomer A to monomer B in an AB co-polymer, shown in \autoref{equation 1} 

\begin{equation}\label{equation 1} 
\mathit{R_{m}}  = \frac{n(A)}{n(A+B)} 
\end{equation} 

Degree of polymerization, DOP, is the mean polymer chain length in the system, shown in \autoref{equation 2}. 

\begin{equation}\label{equation 2} 
DOP = \overline{x}(nA+nB) 
\end{equation} 

Number of chains, n~chains, is the total number of chains built by SwiftPol and appended to the object, shown in \autoref{equation 3}. 

\begin{equation}\label{equation 3} 
n_{chains} = total\,number\,of\,chains\,built 
\end{equation} 

Blockiness, *b*, is a measurement of the distribution of monomers in an AB co-polymer, shown in \autoref{equation 4}. 

\begin{equation}\label{equation 4} 
\mathit{b} = \frac{nB-B\,bonds}{nA-B\,bonds} 
\end{equation} 

Residual monomer, *M~resid*, is the % of residual monomer molecules in the system, shown in equation 5. 

\begin{equation}\label{equation 4} 
\mathit{M_{resid}} = \frac{M_{w}(M_{resid)}{M_{w}(Carbon-containing\,compounds)} 
\end{equation} 
 
# Acknowledgements 
 
Hannah Turney is supported by funding contributions from UKRI Biotechnology and Biological Sciences Research Council (grant ref. BB/T008709/1) and Johnson&Johnson Innovative Medicine. 

 

We acknowledge contributions and feedback from Jeffrey Wagner at the Open Force field consortium and Anusha Lalitha, David Hahn and Gary Tresadern at Johnson&Johnson Innovative Medicine 
 
# References 

Wan, B. et al. (2021) ‘Effect of polymer source on in vitro drug release from PLGA microspheres’, International Journal of Pharmaceutics, 607, p. 120907. https://doi.org/10.1016/j.ijpharm.2021.120907. 

 