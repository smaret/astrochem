"""Python module for Astrochem.

"""

PACKAGE_VERSION = "@PACKAGE_VERSION@"

try:
    import libpyastrochem
    import wrapper
except ImportError:
    pass
import tools
