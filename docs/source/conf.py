# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
# import guzzle_sphinx_theme

sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'gasparse'
copyright = '2024, Francisco Vasconcelos'
author = 'Francisco Vasconcelos'
release = '0.0.3a'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.viewcode',
        'sphinx.ext.napoleon',
        'sphinx.ext.doctest',
        "sphinx.ext.autosectionlabel",
]
autosectionlabel_prefix_document = True


templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


# html_theme_path = guzzle_sphinx_theme.html_theme_path()
# html_permalinks_icon = '<span>#</span>'
# html_theme = 'sphinxawesome_theme'
# html_theme = 'insipid'

# html_theme = 'guzzle_sphinx_theme'
html_theme = 'sphinx_book_theme'

# Register the theme as an extension to generate a sitemap.xml
# extensions.append("guzzle_sphinx_theme")

# Guzzle theme options (see theme.conf for more information)
# html_theme_options = {
    # Set the name of the project to appear in the sidebar
    # "project_nav_name": "gasparse",
    # "globaltoc_collapse": True,
    # "globaltoc_depth": 3,

# }


html_static_path = ['_static']
