# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/stable/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# Incase the project was not installed
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import swiftpol
import os
import sys

# Add the swiftpol module to the Python path
sys.path.insert(0, os.path.abspath('../swiftpol'))

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = 'SwiftPol'
author = 'Hannah Turney'
release = '0.1.3'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'imported-members': False,
    'private-members': False,
    'show-inheritance': True,
}
autodoc_member_order = 'bysource'

# Mock imports for external dependencies
autodoc_mock_imports = [
    'rdkit',
    'numpy',
    'scipy',
    'pandas',
    'matplotlib',
]

# -- Project information -----------------------------------------------------

project = 'swiftpol'
copyright = ("2025, Hannah Turney. Project structure based on the "
             "Computational Molecular Science Python Cookiecutter version 1.1")
author = 'Hannah Turney'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True


templates_path = ['_templates']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'default'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


htmlhelp_basename = 'swiftpoldoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'swiftpol.tex', 'swiftpol Documentation',
     'swiftpol', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'swiftpol', 'swiftpol Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'swiftpol', 'swiftpol Documentation',
     author, 'swiftpol', 'Tools for building polydisperse polymer systems for molecular dynamics.',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------
