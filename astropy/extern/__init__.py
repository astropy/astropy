# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This packages contains python packages that are bundled with Astropy but are
external to Astropy, and hence are developed in a separate source tree.  Note
that this package is distinct from the /cextern directory of the source code
distribution, as that directory only contains C extension code.

See the README.rst in this directory of the Astropy source repository for more
details.
"""

from lazy_loader import attach_stub

__getattr__, __dir__, __all__ = attach_stub(__name__, __file__)
