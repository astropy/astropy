# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from . import FitsTestCase

from ..fitstime import GLOBAL_TIME_INFO, time_to_fits, is_time_column_keyword
from ....coordinates import EarthLocation
from ....table import Table, QTable
from ....time import Time
from ....time.core import BARYCENTRIC_SCALES
from ....time.formats import FITS_DEPRECATED_SCALES
from ....tests.helper import catch_warnings


class TestFitsTime(FitsTestCase):

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
    def test_time_to_fits_loc(self, table_types):
        """
        Test all the unusual conditions for locations of ``Time``
        columns in a ``Table``.
        """
        t = table_types()
        t['a'] = Time(self.time, format='isot', scale='utc')
        t['b'] = Time(self.time, format='isot', scale='tt')

        # Check that vectorized location raises an exception
        t['a'].location = EarthLocation([1,2], [2,3], [3,4])

        with pytest.raises(ValueError) as err:
            table, hdr = time_to_fits(t)
        assert 'Vectorized Location of Time Column' in str(err.value)

        # Check that multiple Time columns with different locations raise an exception
        t['a'].location = EarthLocation(1, 2, 3)
        t['b'].location = EarthLocation(2, 3, 4)

        with pytest.raises(ValueError) as err:
            table, hdr = time_to_fits(t)
            assert 'Multiple Time Columns with different geocentric' in str(err.value)

        # Check that Time column with no location specified will assume global location
        t['b'].location = None

        with catch_warnings() as w:
            table, hdr = time_to_fits(t)
            assert len(w) == 1
            assert str(w[0].message).startswith('Time Column "b" has no specified '
                                                'location, but global Time Position '
                                                'is present')

        # Check that multiple Time columns with same location can be written
        t['b'].location = EarthLocation(1, 2, 3)

        with catch_warnings() as w:
            table, hdr = time_to_fits(t)
            assert len(w) == 0

        # Check compatibility of Time Scales and Reference Positions

        for scale in BARYCENTRIC_SCALES:
            t.replace_column('a', getattr(t['a'], scale))
            with catch_warnings() as w:
                table, hdr = time_to_fits(t)
                assert len(w) == 1
                assert str(w[0].message).startswith('Earth Location "TOPOCENTER" '
                                                    'for Time Column')

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_time_to_fits_header(self, table_types):
        """
        Test the header returned by ``time_to_fits`` explicitly.
        """
        t = table_types()
        t['a'] = Time(self.time, format='isot', scale='utc',
                      location=EarthLocation(-2446354,
                      4237210, 4077985, unit='m'))
        t['b'] = Time([1,2], format='cxcsec', scale='tt')

        ideal_col_hdr = {'OBSGEO-X' : t['a'].location.x.value,
                         'OBSGEO-Y' : t['a'].location.y.value,
                         'OBSGEO-Z' : t['a'].location.z.value}

        table, hdr = time_to_fits(t)

        # Check global time keywords
        for key, value in GLOBAL_TIME_INFO.items():
            assert hdr[key] == value[0]
            assert hdr.comments[key] == value[1]
            hdr.remove(key)

        # Check column-specific(related) keywords
        for key, value in ideal_col_hdr.items():
            assert hdr[key] == value
            hdr.remove(key)

        coord_info = table.meta['__coordinate_columns__']
        for colname in coord_info:
            assert coord_info[colname]['coord_type'] == t[colname].scale.upper()

        assert coord_info['a']['time_ref_pos'] == 'TOPOCENTER'

        assert len(hdr) == 0

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_fits_to_time_meta(self, table_types):
        """
        Test that the relevant global time metadata is read into
        ``Table.meta`` as ``Time``.
        """
        t = table_types()
        t['a'] = Time(self.time, format='isot', scale='utc')
        t.meta['DATE'] = '1999-01-01T00:00:00'
        t.meta['MJD-OBS'] = 56670

        # Test for default read/write behaviour (raw data)
        t.write(self.temp('time.fits'), format='fits', overwrite=True)
        tm = table_types.read(self.temp('time.fits'), format='fits')

        # Test DATE
        assert not isinstance(tm.meta['DATE'], Time)
        assert tm.meta['DATE'] == t.meta['DATE']

        # Test MJD-xxx
        assert not isinstance(tm.meta['MJD-OBS'], Time)
        assert tm.meta['MJD-OBS'] == t.meta['MJD-OBS']

        # Test for native astropy objects read/write behaviour
        t.write(self.temp('time.fits'), format='fits', overwrite=True,
                astropy_native=True)
        tm = table_types.read(self.temp('time.fits'), format='fits',
                              astropy_native=True)

        # Test DATE
        assert isinstance(tm.meta['DATE'], Time)
        assert tm.meta['DATE'].value == t.meta['DATE'] + '(UTC)'
        assert tm.meta['DATE'].format == 'fits'
        # Default time scale according to the FITS standard is UTC
        assert tm.meta['DATE'].scale == 'utc'

        # Test MJD-xxx
        assert isinstance(tm.meta['MJD-OBS'], Time)
        assert tm.meta['MJD-OBS'].value == t.meta['MJD-OBS']
        assert tm.meta['MJD-OBS'].format == 'mjd'
        assert tm.meta['MJD-OBS'].scale == 'utc'

        # Explicitly specified Time Scale
        t.meta['TIMESYS'] = 'ET'

        t.write(self.temp('time.fits'), format='fits', overwrite=True,
                astropy_native=True)
        tm = table_types.read(self.temp('time.fits'), format='fits',
                              astropy_native=True)

        # Test DATE
        assert isinstance(tm.meta['DATE'], Time)
        assert tm.meta['DATE'].value == t.meta['DATE'] + '(UTC)'
        assert tm.meta['DATE'].scale == 'utc'

        # Test MJD-xxx
        assert isinstance(tm.meta['MJD-OBS'], Time)
        assert tm.meta['MJD-OBS'].value == t.meta['MJD-OBS']
        assert tm.meta['MJD-OBS'].scale == FITS_DEPRECATED_SCALES[t.meta['TIMESYS']]
