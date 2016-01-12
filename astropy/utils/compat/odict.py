# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import

import warnings
from collections import OrderedDict as _OrderedDict

from ..exceptions import AstropyDeprecationWarning


class OrderedDict(_OrderedDict):
    def __init__(self, *args, **kwargs):
        warnings.warn("astropy.utils.compat.odict.OrderedDict is now deprecated - import OrderedDict from the collections module instead", AstropyDeprecationWarning)
        super(OrderedDict, self).__init__(*args, **kwargs)
