# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
# Project information

project = "Fractal Client"
copyright = (
    "2022, Friedrich Miescher Institute for Biomedical Research and "
    "University of Zurich"
)
version = "1.3.0a1"
language = "en"

# General configuration

extensions = [
    "sphinxarg.ext",
    "sphinx_rtd_theme",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinxarg.ext",
]

# Options for HTML output

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "logo_only": False,
    "sticky_navigation": True,
    "titles_only": True,
    "navigation_depth": 5,
    "prev_next_buttons_location": None,
    "display_version": True,
    "style_external_links": True,
}

# Disable python default syntax highlighting (see
# https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-highlight_language)
highlight_language = "none"
