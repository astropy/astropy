# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import

import warnings
from collections import OrderedDict

from ..exceptions import AstropyDeprecationWarning

warnings.warn("astropy.utils.compat.odict is now deprecated - import OrderedDict from the collections module instead")

