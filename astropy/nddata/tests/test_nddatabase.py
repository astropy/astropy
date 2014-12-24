# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Tests of NDDataBase

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..nddatabase import NDDataBase
from ...tests.helper import pytest


class MinimalSubclass(NDDataBase):
    def __init__(self):
        super(MinimalSubclass, self).__init__()

    @property
    def data(self):
        return None

    @property
    def mask(self):
        return super(MinimalSubclass, self).mask

    @property
    def unit(self):
        return super(MinimalSubclass, self).unit

    @property
    def wcs(self):
        return super(MinimalSubclass, self).wcs

    @property
    def meta(self):
        return super(MinimalSubclass, self).meta


class MinimalUncertainty(object):
    """
    Define the minimum attributes acceptable as an uncertainty object.
    """
    def __init__(self, value):
        self._uncertainty = value

    @property
    def uncertainty_type(self):
        return "totally and completely fake"


def test_nddatabase_subclass():
    a = MinimalSubclass()
    assert a.meta is None
    assert a.data is None
    assert a.mask is None
    assert a.unit is None
    assert a.wcs is None
    good_uncertainty = MinimalUncertainty(5)
    a.uncertainty = good_uncertainty
    assert a.uncertainty is good_uncertainty
    bad_uncertainty = 5
    with pytest.raises(TypeError):
        a.uncertainty = bad_uncertainty
