project = "yeti-iga"

copyright = '2025, LaMCoS'
author = 'Arnaud Duval'


try:
    from importlib.metadata import version as get_version
except ImportError:
    from importlib_metadata import version as get_version  # Python < 3.8

release = get_version(project)
version = ".".join(release.split(".")[:2])  # short version, ex: "1.2"

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
    'version_selector': True,
}


bibtex_bibfiles = [
    'YETI.bib',
]