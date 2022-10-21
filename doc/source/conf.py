import sphinx_rtd_theme

project = "yeti"

copyright = '@SPHINX_TARGET_YEAR@, LaMCoS'
author = 'Arnaud Duval'

version = '@SPHINX_TARGET_VERSION_MAJOR@.@SPHINX_TARGET_VERSION_MINOR@'
release = '@SPHINX_TARGET_VERSION@'


extensions = [
    'sphinx_rtd_theme',
    'sphinxcontrib.bibtex',
]

templates_path = ['_templates']

exclude_patterns = []
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

bibtex_bibfiles = [
    'YETI.bib',
]