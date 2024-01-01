import os
import sys
current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../../"))
sys.path.insert(0, target_dir)

#import mock
 
#MOCK_MODULES = ['numpy', 'scipy', 'mpmath']
#for mod_name in MOCK_MODULES:
#	sys.modules[mod_name] = mock.Mock()

print(target_dir)

# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Algebraic differentiators'
copyright = '2023, Amine Othmane'
author = 'Amine Othmane'

# The full version, including alpha/beta/rc tags
release = '2.1'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme'
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    "body_max_width": "none"
}

# -- Options for EPUB output
epub_show_urls = 'footnote'
