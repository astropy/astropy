# coding: utf-8
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.testing import assert_array_equal
import pytest

from astropy import units as u
from astropy.table import QTable

from astropy.utils.masked import Masked
from astropy.utils.compat.optional_deps import HAS_YAML, HAS_H5PY

from .test_masked import assert_masked_equal


FILE_FORMATS = ['ecsv', 'fits']
if HAS_H5PY:
    FILE_FORMATS.append('h5')


class MaskedArrayTableSetup:
    @classmethod
    def setup_arrays(self):
        self.a = np.array([3., 5., 0.])
        self.mask_a = np.array([True, False, False])

    @classmethod
    def setup_class(self):
        self.setup_arrays()
        self.ma = Masked(self.a, mask=self.mask_a)
        self.ma.info.format = '.1f'
        self.t = QTable([self.ma], names=['ma'])


class TestMaskedArrayTable(MaskedArrayTableSetup):
    def test_table_initialization(self):
        assert_array_equal(self.t['ma'].unmasked, self.a)
        assert_array_equal(self.t['ma'].mask, self.mask_a)
        assert repr(self.t).splitlines()[-3:] == [
            '    ———',
            '    5.0',
            '    0.0']

    @pytest.mark.skipif(not HAS_YAML, reason='serialization needs yaml')
    @pytest.mark.parametrize('file_format', FILE_FORMATS)
    def test_table_write(self, file_format, tmpdir):
        name = str(tmpdir.join(f"a.{file_format}"))
        kwargs = {}
        if file_format == 'h5':
            kwargs['path'] = 'trial'
            kwargs['serialize_meta'] = True

        self.t.write(name, **kwargs)
        t2 = QTable.read(name)
        assert isinstance(t2['ma'], self.ma.__class__)
        assert np.all(t2['ma'] == self.ma)
        assert np.all(t2['ma'].mask == self.mask_a)
        if file_format == 'fits' and type(self.a) is np.ndarray:
            # Imperfect roundtrip through FITS native format description.
            assert self.t['ma'].info.format in t2['ma'].info.format
        else:
            assert t2['ma'].info.format == self.t['ma'].info.format


@pytest.mark.skipif(not HAS_YAML, reason='serialization needs yaml')
class TestSerializationMethods(MaskedArrayTableSetup):
    # TODO: ensure this works for MaskedQuantity, etc., as well.
    # Needs to somehow pass on serialize_method; see MaskedArraySubclassInfo.
    @pytest.mark.parametrize('serialize_method', ['data_mask', 'null_value'])
    def test_table_write(self, serialize_method, tmpdir):
        name = str(tmpdir.join("test.ecsv"))
        self.t.write(name, serialize_method=serialize_method)
        with open(name) as fh:
            lines = fh.readlines()

        t2 = QTable.read(name)
        assert isinstance(t2['ma'], self.ma.__class__)

        if serialize_method == 'data_mask':
            # Will data_mask, we have data and mask columns and should
            # exactly round-trip.
            assert len(lines[-1].split()) == 2
            assert_masked_equal(t2['ma'], self.ma)
        else:
            # With null_value we have just a data column with null values
            # marked, so we lost information on the data below the mask.
            assert len(lines[-1].split()) == 1
            assert np.all(t2['ma'] == self.ma)
            assert np.all(t2['ma'].mask == self.mask_a)

    def test_non_existing_serialize_method(self, tmpdir):
        name = str(tmpdir.join('bad.ecsv'))
        with pytest.raises(ValueError, match='serialize method must be'):
            self.t.write(name, serialize_method='bad_serialize_method')


class TestMaskedQuantityTable(TestMaskedArrayTable):
    @classmethod
    def setup_arrays(self):
        self.a = np.array([3., 5., 0.]) << u.m
        self.mask_a = np.array([True, False, False])
