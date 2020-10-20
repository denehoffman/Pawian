"""Configuration file for the Sphinx documentation builder.

This file only contains a selection of the most common options. For a full
list see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

# -- Project information -----------------------------------------------------
project = "PAWIAN"
copyright = "2020, Ruhr-Universtät Bochum"
author = "Bertram Kopf"


# -- General configuration ---------------------------------------------------
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

# The master toctree document.
master_doc = "index"
rst_prolog = open("global.inc", "r").read()

extensions = [
    "myst_parser",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.ifconfig",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
]
exclude_patterns = [
    "**.ipynb_checkpoints",
    "*build",
]

# General sphinx settings
add_module_names = False
highlight_language = "shell"
html_copy_source = False  # do not copy rst files
html_show_copyright = False
html_show_sourcelink = False
html_show_sphinx = False
html_theme = "sphinx_rtd_theme"
pygments_style = "sphinx"
todo_include_todos = False
viewcode_follow_imported_members = True

# Cross-referencing configuration
# default_role = "cpp"
# primary_domain = "cpp"
nitpicky = True  # warn if cross-references are missing
nitpick_ignore = []

# Intersphinx settings
intersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
    "python": ("https://docs.python.org/3", None),
}

# Settings for autosectionlabel
autosectionlabel_prefix_document = True

# Settings for linkcheck
linkcheck_anchors = False

# Settings for myst-parser
myst_update_mathjax = False
