# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from . import FitsTestCase

from ..fitstime import GLOBAL_TIME_INFO, is_time_column_keyword, Time_to_FITS
from ....coordinates import EarthLocation
from ....table import Table, QTable
from ....time import Time, BARYCENTRIC_SCALES
from ....tests.helper import catch_warnings


class Test_FITStime(FitsTestCase):

    def setup_class(self):
        self.time = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00']

    def test_is_time_column_keyword(self):
        # Time column keyword without column number
        assert is_time_column_keyword('TRPOS') is False

        # Global time column keyword
        assert is_time_column_keyword('TIMESYS') is False

        # Valid time column keyword
        assert is_time_column_keyword('TRPOS12') is True

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_Time_to_FITS_loc(self, table_types):
        """
        Test all the unusual conditions for locations of `Time` columns in a `Table`.
        """
        t = table_types()
        t['a'] = Time(self.time, format='isot', scale='utc')
        t['b'] = Time(self.time, format='isot', scale='tt')

        # Check that vectorized location raises an exception
        t['a'].location = EarthLocation([1,2], [2,3], [3,4,])

        with pytest.raises(ValueError) as err:
            table, hdr = Time_to_FITS(t)
        assert 'Vectorized Location of Time Column' in str(err.value)

        # Check that multiple Time columns with different locations raise an exception
        t['a'].location = EarthLocation(1, 2, 3)
        t['b'].location = EarthLocation(2, 3, 4)

        with pytest.raises(ValueError) as err:
            table, hdr = Time_to_FITS(t)
            assert 'Multiple Time Columns with different geocentric' in str(err.value)

        # Check that Time column with no location specified will assume global location
        t['b'].location = None

        with catch_warnings() as w:
            table, hdr = Time_to_FITS(t)
            assert len(w) == 1
            assert str(w[0].message).startswith('Time Column "b" has no specified location,'
                                                ' but global Time Position is present')

        # Check that multiple Time columns with same location can be written
        t['b'].location = EarthLocation(1, 2, 3)

        with catch_warnings() as w:
            table, hdr = Time_to_FITS(t)
            assert len(w) == 0

        # Check compatibility of Time Scales and Reference Positions

        for scale in BARYCENTRIC_SCALES:
            t.replace_column('a', getattr(t['a'], scale))
            with catch_warnings() as w:
                table, hdr = Time_to_FITS(t)
                assert len(w) == 1
                assert str(w[0].message).startswith('Earth Location "TOPOCENTER" for Time Column')

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_Time_to_FITS_header(self, table_types):
        """
        Test the header returned by `Time_to_FITS` explicitly.
        """
        t = table_types()
        t['a'] = Time(self.time, format='isot', scale='utc',
                      location=EarthLocation(-2446353.80003635,
                      4237209.07495215, 4077985.57220038, unit='m'))
        t['b'] = Time([1,2], format='cxcsec', scale='tt')

        ideal_col_hdr = {'OBSGEO-X' : t['a'].location.x.value,
                         'OBSGEO-Y' : t['a'].location.y.value,
                         'OBSGEO-Z' : t['a'].location.z.value,
                         'TCTYP1' : t['a'].scale.upper(),
                         'TRPOS1' : 'TOPOCENTER',
                         'TCTYP2' : t['b'].scale.upper()}

        table, hdr = Time_to_FITS(t)

        # Check global time keywords
        for key, value in GLOBAL_TIME_INFO.items():
            assert hdr[key] == value[0]
            assert hdr.comments[key] == value[1]
            hdr.remove(key)

        # Check column-specific(related) keywords
        for key, value in ideal_col_hdr.items():
            assert hdr[key] == value
            hdr.remove(key)

        assert len(hdr) == 0

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_FITS_to_Time_meta(self, table_types):
        """
        Test that the relevant global time metadata is read into `Table.meta` as `Time`.
        """
        t = Table()
        t['a'] = Time(self.time, format='isot', scale='utc')
        t.meta['DATE'] = '1999-01-01T00:00:00'

        t.write(self.temp('time.fits'), format='fits')

        tm = Table.read(self.temp('time.fits'), format='fits')

        assert isinstance(tm.meta['DATE'], Time)
        assert tm.meta['DATE'].value == t.meta['DATE']
        # Default time scale according to the FITS standard is UTC
        assert tm.meta['DATE'].scale == 'utc'

        # Explicitly specified Time Scale
        t.meta['TIMESYS'] = 'TT'

        t.write(self.temp('time.fits'), format='fits', overwrite=True)

        tm = Table.read(self.temp('time.fits'), format='fits')

        assert isinstance(tm.meta['DATE'], Time)
        assert tm.meta['DATE'].value == t.meta['DATE']
        assert tm.meta['DATE'].scale == t.meta['TIMESYS'].lower()
