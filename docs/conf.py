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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import subprocess, os

#read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

#if read_the_docs_build:
subprocess.call('cd ../; doxygen docs/Doxyfile', shell=True)

# -- Project information -----------------------------------------------------

project = 'QSL'
copyright = '2020, Lana Mineh and John Scott'
author = 'Lana Mineh and John Scott'

# Replace |project| with the project name above in docs
rst_epilog = '.. |project| replace:: %s' % project

# The full version, including alpha/beta/rc tags
release = '0.1'

# This will print the version under the program title
#version = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'breathe', 'sphinx.ext.todo', 'sphinx.ext.mathjax', 'sphinx_copybutton',
]

# Breathe configuration
breathe_projects = { "QSL": "doxy/xml/", }
breathe_default_project = "QSL"
todo_include_todos = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
#html_theme = 'alabaster'

copybutton_prompt_text = "$ "

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []
