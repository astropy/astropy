# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This subpackage contains utilities used by Astropy. Public functions and class
here are same to be used by affiliated or other packages, but this package is
for utilities that are mostly of use to developers.  This also includes the
`astropy.utils.compat` package, which houses utilities that provide 
compatibility and bugfixes across all versions of Python that Astropy supports.

For astronomy-related "utilities" of general use (e.g. not specific to some 
other subpackage), see the `astropy.tools` package. 
"""

from .compat.odict import OrderedDict
