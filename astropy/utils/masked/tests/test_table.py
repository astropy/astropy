# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from astropy import units as u
from astropy.table import QTable

from .. import Masked

try:
    import yaml  # noqa
except ImportError:
    FILE_FORMATS = []
else:
    # FITS cannot yet deal with float masked columns generally.
    FILE_FORMATS = ['ecsv']
    try:
        import h5py  # noqa
        FILE_FORMATS.append('h5')
    except ImportError:
        pass


class TestMaskedArrayTable:
    @classmethod
    def setup_arrays(self):
        self.a = np.array([3., 5., 0.])
        self.mask_a = np.array([True, False, False])

    def setup_class(self):
        self.setup_arrays()
        self.ma = Masked(self.a, mask=self.mask_a)
        self.ma.info.format = '.1f'
        self.t = QTable([self.ma], names=['ma'])

    def test_table_initialization(self):
        assert_array_equal(self.t['ma'].unmasked, self.a)
        assert_array_equal(self.t['ma'].mask, self.mask_a)
        assert repr(self.t).splitlines()[-3:] == [
            '    ———',
            '    5.0',
            '    0.0']

    @pytest.mark.parametrize('file_format', FILE_FORMATS)
    def test_table_write(self, file_format, tmpdir):
        name = str(tmpdir / f"a.{file_format}")
        if file_format == 'h5':
            self.t.write(name, path='trial', serialize_meta=True)
        else:
            self.t.write(name)
        t2 = QTable.read(name)
        assert isinstance(t2['ma'], self.ma.__class__)
        assert np.all(t2['ma'] == self.ma)
        assert np.all(t2['ma'].mask == self.mask_a)


class TestMaskedQuantityTable(TestMaskedArrayTable):
    @classmethod
    def setup_arrays(self):
        self.a = np.array([3., 5., 0.]) << u.m
        self.mask_a = np.array([True, False, False])
