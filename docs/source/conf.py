import tomli
import warnings


# Set version from pyproject.toml file
with open('../../pyproject.toml', 'rb') as f:
    toml_data = tomli.load(f)
    release = toml_data['project']['version']

project = "yeti"

copyright = '2024, LaMCoS'
author = 'Arnaud Duval'

# version = '@SPHINX_TARGET_VERSION_MAJOR@.@SPHINX_TARGET_VERSION_MINOR@'
# release = '@SPHINX_TARGET_VERSION@'
version = release

numfig = True

extensions = [
    'sphinx.ext.autosectionlabel',
    'sphinx_rtd_theme',
    'sphinxcontrib.bibtex',
    'sphinx.ext.mathjax',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon'       # Handle Numpy or Google style docstrings
]

templates_path = ['_templates']

exclude_patterns = []
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/logo-yeti_with_text.svg'
html_theme_option = {
    'logo_only': True,
    'display_version': False,
}

bibtex_bibfiles = [
    'YETI.bib',
]