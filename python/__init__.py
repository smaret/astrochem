"""Python module for Astrochem.

"""

PACKAGE_VERSION = "@PACKAGE_VERSION@"

try:
    from . import libpyastrochem
    from . import wrapper
except ImportError:
    pass
from . import tools
