# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Tests of NDDataBase

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ..nddata_base import NDDataBase


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

    @property
    def uncertainty(self):
        return super(MinimalSubclass, self).uncertainty


def test_nddata_base_subclass():
    a = MinimalSubclass()
    assert a.meta is None
    assert a.data is None
    assert a.mask is None
    assert a.unit is None
    assert a.wcs is None
    assert a.uncertainty is None
