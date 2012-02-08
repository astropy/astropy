# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage contains developer-oriented utilities used by Astropy.

Public functions and classes in this subpackage are safe to be used by other
packages, but this subpackage is for utilities that are primarily of use for
developers or to implement python hacks. This subpackage also includes the
`astropy.utils.compat` package, which houses utilities that provide
compatibility and bugfixes across all versions of Python that Astropy supports.

For astronomy-specific utilities of general use (e.g. not specific to some
other subpackage), see the `astropy.tools` package.
"""

from .compat.odict import OrderedDict
from .misc import *

# The location of the online documentation for astropy
# TODO: For now, this is just the latest version at readthedocs.org
# In the future this should correspond to the current version
# of astropy.
online_docs_root = 'http://astropy.readthedocs.org/en/latest/'
