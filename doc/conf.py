# Configuration file for the Sphinx documentation builder.

import sys
import os
from datetime import datetime

sys.path.insert(0, os.path.abspath("../python/"))

# Mock out the import of Numpy, to compile on ReadTheDoc
from unittest.mock import Mock as MagicMock


class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return Mock()


MOCK_MODULES = ["numpy", "libpyastrochem"]
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

import tools  # noqa: E402, F401
import wrapper  # noqa: E402, F401

project = "Astrochem"
author = "Sébastien Maret"
thisyear = datetime.today().year
copyright = f"2006-{thisyear}, {author}"
release = "0.10"

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
]

html_theme = "sphinx_rtd_theme"
