# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
#import dolfen

sys.path.insert(0, os.path.abspath("../dolfen/"))

# -- Project information -----------------------------------------------------

project = "dolfen"
copyright = "2023, Jethro Linley"
author = "Jethro Linley"

# The full version, including alpha/beta/rc tags
#release = dolfen.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    "sphinx.ext.inheritance_diagram",
    "autoapi.extension",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# -- Configure autoapi -------------------------------------------------------
autoapi_type = "python"
autoapi_dirs = ["../dolfen/"]
autoapi_add_toctree_entry = False
autoapi_options = ["members", "show-inheritance", "show-module-summary"]

# -- RST config --------------------------------------------------------------
# Inline python highlighting, base on https://stackoverflow.com/a/68931907
rst_prolog = """
.. role:: python(code)
    :language: python
    :class: highlight
"""
