from pathlib import Path
from typing import Any
from typing import List

from sphinx.application import Sphinx
from sphinx.ext import apidoc

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Fractal Client'
copyright = (
    "2022, Friedrich Miescher Institute for Biomedical Research and "
    "University of Zurich"
)
version = "1.0.0"
language = "en"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinxarg.ext',
    "sphinx_rtd_theme",
]

templates_path = ['_templates']
exclude_patterns = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinxarg.ext',

]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    "logo_only": False,
    "sticky_navigation": True,
    "titles_only": True,
    "navigation_depth": 5,
    "prev_next_buttons_location": None,
    "display_version": True,
    "style_external_links": True,
}
