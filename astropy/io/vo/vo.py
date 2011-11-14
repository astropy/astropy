import warnings
warnings.warn(
    "vo is deprecated.  Use astropy.io.vo instead.",
    DeprecationWarning)

from astropy.io.vo import *
from astropy.io.vo import converters, table, tree, ucd, util, xmlutil
from astropy.io.vo import exceptions as voexceptions
