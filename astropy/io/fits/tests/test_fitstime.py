# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

import pytest
import numpy as np

from . import FitsTestCase

from ..fitstime import GLOBAL_TIME_INFO, time_to_fits, is_time_column_keyword
from ....coordinates import EarthLocation
from ....io import fits
from ....table import Table, QTable
from ....time import Time, TimeDelta
from ....time.core import BARYCENTRIC_SCALES
from ....time.formats import FITS_DEPRECATED_SCALES
from ....tests.helper import catch_warnings


class TestFitsTime(FitsTestCase):

    def setup_class(self):
        self.time = np.array(['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00'])

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
        Test the header and metadata returned by ``time_to_fits``.
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

        # Check the global time keywords in hdr
        for key, value in GLOBAL_TIME_INFO.items():
            assert hdr[key] == value[0]
            assert hdr.comments[key] == value[1]
            hdr.remove(key)

        for key, value in ideal_col_hdr.items():
            assert hdr[key] == value
            hdr.remove(key)

        # Check the column-specific time metadata
        coord_info = table.meta['__coordinate_columns__']
        for colname in coord_info:
            assert coord_info[colname]['coord_type'] == t[colname].scale.upper()
            assert coord_info[colname]['coord_unit'] == 'd'

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

        # Test for default write behaviour (full precision) and read it
        # back using native astropy objects; thus, ensure its round-trip
        t.write(self.temp('time.fits'), format='fits', overwrite=True)
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

        t.write(self.temp('time.fits'), format='fits', overwrite=True)
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

        # Test for conversion of time data to its value, as defined by its format
        t['a'].info.serialize_method['fits'] = 'formatted_value'
        t.write(self.temp('time.fits'), format='fits', overwrite=True)
        tm = table_types.read(self.temp('time.fits'), format='fits')

        # Test DATE
        assert not isinstance(tm.meta['DATE'], Time)
        assert tm.meta['DATE'] == t.meta['DATE']

        # Test MJD-xxx
        assert not isinstance(tm.meta['MJD-OBS'], Time)
        assert tm.meta['MJD-OBS'] == t.meta['MJD-OBS']

        assert (tm['a'] == t['a'].value).all()

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_time_loc_unit(self, table_types):
        """
        Test that ``location`` specified by using any valid unit
        (length/angle) in ``Time`` columns gets stored in FITS
        as ITRS Cartesian coordinates (X, Y, Z), each in m.
        Test that it round-trips through FITS.
        """
        t = table_types()
        t['a'] = Time(self.time, format='isot', scale='utc',
                      location=EarthLocation(1,2,3, unit='km'))

        table, hdr = time_to_fits(t)

        # Check the header
        hdr['OBSGEO-X'] == t['a'].location.x.to_value(unit='m')
        hdr['OBSGEO-Y'] == t['a'].location.y.to_value(unit='m')
        hdr['OBSGEO-Z'] == t['a'].location.z.to_value(unit='m')

        t.write(self.temp('time.fits'), format='fits', overwrite=True)
        tm = table_types.read(self.temp('time.fits'), format='fits',
                              astropy_native=True)

        # Check the round-trip of location
        tm['a'].location == t['a'].location
        tm['a'].location.x.value == t['a'].location.x.to_value(unit='m')
        tm['a'].location.y.value == t['a'].location.y.to_value(unit='m')
        tm['a'].location.z.value == t['a'].location.z.to_value(unit='m')

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_io_time_read_fits(self, table_types):
        """
        Test that FITS table with time columns (standard compliant)
        can be read by io.fits as a table with Time columns.
        This tests the following:
        1. The special-case where a column has the name 'TIME' and a
           time unit
        2. Time from Epoch (Reference time) is appropriately converted.
        3. Coordinate columns (corresponding to coordinate keywords in the header)
           other than time, that is, spatial coordinates, are not mistaken
           to be time.
        """
        filename = self.data('chandra_time.fits')
        tm = table_types.read(filename, astropy_native=True)

        # Test case 1
        assert isinstance(tm['time'], Time)
        assert tm['time'].scale == 'tt'
        assert tm['time'].format == 'mjd'

        non_native = table_types.read(filename)

        # Test case 2
        ref_time = Time(non_native.meta['MJDREF'], format='mjd',
                        scale=non_native.meta['TIMESYS'].lower())
        delta_time = TimeDelta(non_native['time'])
        assert (ref_time + delta_time == tm['time']).all()

        # Test case 3
        for colname in ['chipx', 'chipy', 'detx', 'dety', 'x', 'y']:
            assert not isinstance(tm[colname], Time)

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_io_time_read_fits_datetime(self, table_types):
        """
        Test that ISO-8601 Datetime String Columns are read correctly.
        """
        # Datetime column
        c = fits.Column(name='datetime', format='A29', coord_type='TCG',
                        time_ref_pos='GEOCENTER', array=self.time)

        # Explicitly create a FITS Binary Table
        bhdu = fits.BinTableHDU.from_columns([c])
        bhdu.writeto(self.temp('time.fits'), overwrite=True)

        tm = table_types.read(self.temp('time.fits'), astropy_native=True)

        assert isinstance(tm['datetime'], Time)
        assert tm['datetime'].scale == 'tcg'
        assert tm['datetime'].format == 'fits'
        assert (tm['datetime'] == self.time).all()

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_io_time_read_fits_location(self, table_types):
        """
        Test that geocentric/geodetic observatory position is read
        properly, as and when it is specified.
        """
        # Datetime column
        c = fits.Column(name='datetime', format='A29', coord_type='TT',
                        time_ref_pos='TOPOCENTER', array=self.time)

        # Observatory position in ITRS Cartesian coordinates (geocentric)
        cards = [('OBSGEO-X', -2446354), ('OBSGEO-Y', 4237210),
                 ('OBSGEO-Z', 4077985)]

        # Explicitly create a FITS Binary Table
        bhdu = fits.BinTableHDU.from_columns([c], header=fits.Header(cards))
        bhdu.writeto(self.temp('time.fits'), overwrite=True)

        tm = table_types.read(self.temp('time.fits'), astropy_native=True)

        assert isinstance(tm['datetime'], Time)
        assert tm['datetime'].location.x.value == -2446354
        assert tm['datetime'].location.y.value == 4237210
        assert tm['datetime'].location.z.value == 4077985

        # Observatory position in geodetic coordinates
        cards = [('OBSGEO-L', 0), ('OBSGEO-B', 0), ('OBSGEO-H', 0)]

        # Explicitly create a FITS Binary Table
        bhdu = fits.BinTableHDU.from_columns([c], header=fits.Header(cards))
        bhdu.writeto(self.temp('time.fits'), overwrite=True)

        tm = table_types.read(self.temp('time.fits'), astropy_native=True)

        assert isinstance(tm['datetime'], Time)
        assert tm['datetime'].location.lon.value == 0
        assert tm['datetime'].location.lat.value == 0
        assert np.isclose(tm['datetime'].location.height.value, 0,
                          rtol=0, atol=1e-9)

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_io_time_read_fits_scale(self, table_types):
        """
        Test handling of 'GPS' and 'LOCAL' time scales which are
        recognized by the FITS standard but are not native to astropy.
        """
        # GPS scale column
        gps_time = np.array([630720013, 630720014])
        c = fits.Column(name='gps_time', format='D', unit='s', coord_type='GPS',
                        coord_unit='s', time_ref_pos='TOPOCENTER', array=gps_time)

        cards = [('OBSGEO-L', 0), ('OBSGEO-B', 0), ('OBSGEO-H', 0)]

        bhdu = fits.BinTableHDU.from_columns([c], header=fits.Header(cards))
        bhdu.writeto(self.temp('time.fits'), overwrite=True)

        with catch_warnings() as w:
            tm = table_types.read(self.temp('time.fits'), astropy_native=True)
            assert len(w) == 1
            assert 'FITS recognized time scale value "GPS"' in str(w[0].message)

        assert isinstance(tm['gps_time'], Time)
        assert tm['gps_time'].format == 'gps'
        assert tm['gps_time'].scale == 'tai'
        assert (tm['gps_time'].value == gps_time).all()

        # LOCAL scale column
        local_time = np.array([1, 2])
        c = fits.Column(name='local_time', format='D', unit='d',
                        coord_type='LOCAL', coord_unit='d',
                        time_ref_pos='RELOCATABLE', array=local_time)

        bhdu = fits.BinTableHDU.from_columns([c])
        bhdu.writeto(self.temp('time.fits'), overwrite=True)

        with catch_warnings() as w:
            tm = table_types.read(self.temp('time.fits'), astropy_native=True)
            assert len(w) == 1
            assert 'FITS recognized time scale value "LOCAL"' in str(w[0].message)

        assert isinstance(tm['local_time'], Time)
        assert tm['local_time'].format == 'mjd'
        # Default scale is UTC
        assert tm['local_time'].scale == 'utc'
        assert (tm['local_time'].value == local_time).all()

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_io_time_read_fits_location_warnings(self, table_types):
        """
        Test warnings for time column reference position.
        """
        # Time reference position "TOPOCENTER" without corresponding
        # observatory position.
        c = fits.Column(name='datetime', format='A29', coord_type='TT',
                        time_ref_pos='TOPOCENTER', array=self.time)

        bhdu = fits.BinTableHDU.from_columns([c])
        bhdu.writeto(self.temp('time.fits'), overwrite=True)

        with catch_warnings() as w:
            tm = table_types.read(self.temp('time.fits'), astropy_native=True)
            assert len(w) == 1
            assert ('observatory position is not properly specified' in
                    str(w[0].message))

        # Default value for time reference position is "TOPOCENTER"
        c = fits.Column(name='datetime', format='A29', coord_type='TT',
                        array=self.time)

        bhdu = fits.BinTableHDU.from_columns([c])
        bhdu.writeto(self.temp('time.fits'), overwrite=True)
        with catch_warnings() as w:
            tm = table_types.read(self.temp('time.fits'), astropy_native=True)
            assert len(w) == 1
            assert ('"TRPOSn" is not specified. The default value for '
                    'it is "TOPOCENTER"' in str(w[0].message))
