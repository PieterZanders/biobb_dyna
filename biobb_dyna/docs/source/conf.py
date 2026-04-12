# -*- coding: utf-8 -*-
# Sphinx configuration for biobb_dyna (aligned with biobb_pytorch docs layout).

import sys
from pathlib import Path

sys.path.insert(0, str(Path("../../").resolve()))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "recommonmark",
]

napoleon_numpy_docstring = False
napoleon_google_docstring = True

templates_path = ["_templates"]
source_suffix = [".rst", ".md"]
master_doc = "index"

project = "biobb_dyna"
copyright = "2026, Bioexcel Project"
author = "Bioexcel Project"

version = "5.0.0"
release = "5.0.0"

language = "en"

exclude_patterns: list[str] = []

pygments_style = "sphinx"
todo_include_todos = False


def setup(app):
    app.add_css_file("theme_overrides.css")
    app.add_js_file("theme_overrides.js")


html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
htmlhelp_basename = "biobb_dyna_doc"

latex_elements: dict[str, str] = {}

latex_documents = [
    (
        master_doc,
        "biobb_dyna.tex",
        "biobb_dyna Documentation",
        "Bioexcel Project",
        "manual",
    ),
]

man_pages = [(master_doc, "biobb_dyna", "biobb_dyna Documentation", [author], 1)]

texinfo_documents = [
    (
        master_doc,
        "biobb_dyna",
        "biobb_dyna Documentation",
        author,
        "biobb_dyna",
        "biobb_dyna MD correlation and network building blocks",
        "Miscellaneous",
    ),
]
