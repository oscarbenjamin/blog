# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'blog'
copyright = '2023, Oscar Benjamin'
author = 'Oscar Benjamin'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

extensions = [
    'sphinx.ext.doctest',
    'sphinx.ext.graphviz',
]
graphviz_output_format = 'svg'

# https://sphinx-comments.readthedocs.io/en/latest/utterances.html
comments_config = {
   "utterances": {
      "repo": "oscarbenjamin/blog",
      "optional": "config",
   }
}
