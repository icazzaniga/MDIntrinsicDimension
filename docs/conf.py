# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from pathlib import Path

import tomllib


def _read_version() -> str:
    """Read the package version from pyproject.toml for docs coherency."""
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    try:
        with pyproject_path.open("rb") as handle:
            return tomllib.load(handle)["project"]["version"]
    except Exception:
        return "0.0.0"


project = 'MDIntrinsicDimension'
copyright = '2025, Irene Cazzaniga'
author = 'Irene Cazzaniga'
release = _read_version()
version = release

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_rtd_theme',
    'myst_nb',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'English'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


html_theme = 'sphinx_rtd_theme'
html_static_path = ["_static"]
html_theme_options = {
    'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#ceff99',
    'flyout_display': 'hidden',
    'version_selector': True,
    'language_selector': True,
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': False,
    'navigation_depth': 2,
    'includehidden': True,
    'titles_only': False,
}
 
html_static_path = ["_static"]
html_css_files = ["custom.css"]


nb_execution_mode = 'off'

nb_remove_input_tags = ["remove-input"]
nb_remove_output_tags = ["remove-output"]
nb_remove_cell_tags = ["remove-cell"]
    
rst_prolog = """
.. role:: mark(code)
   :class: mark
"""

exclude_patterns = ["examples", "maintainer"]
