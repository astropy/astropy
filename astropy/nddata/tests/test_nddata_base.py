# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Tests of NDDataBase


from astropy.nddata.nddata_base import NDDataBase


class MinimalSubclass(NDDataBase):
    def __init__(self):
        super().__init__()

    @property
    def data(self):
        return None

    @property
    def mask(self):
        return super().mask

    @property
    def unit(self):
        return super().unit

    @property
    def wcs(self):
        return super().wcs

    @property
    def meta(self):
        return super().meta

    @property
    def uncertainty(self):
        return super().uncertainty


def test_nddata_base_subclass():
    a = MinimalSubclass()
    assert a.meta is None
    assert a.data is None
    assert a.mask is None
    assert a.unit is None
    assert a.wcs is None
    assert a.uncertainty is None
