#This is the configuration file for the pyfits namespace.

from __future__ import division # confidence high

try:
    import pkg_resources
    __version__ = pkg_resources.get_distribution('pyfits').version
except:
    __version__ = ''

# Import the pyfits core module.
import pyfits.core
import pyfits.util
from pyfits.core import *
from pyfits.util import *

__doc__ = pyfits.core.__doc__

__all__ = pyfits.core.__all__ + pyfits.util.__all__

try:
    import stsci.tools.tester
    def test(*args,**kwds):
        stsci.tools.tester.test(modname=__name__, *args, **kwds)
except ImportError:
    pass

