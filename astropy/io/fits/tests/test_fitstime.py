# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from . import FitsTestCase

from ..fitstime import FITS_time
from ....coordinates import EarthLocation
from ....table import Table, QTable
from ....time import Time
from ....tests.helper import catch_warnings


class Test_FITS_time(FitsTestCase):

    def setup_class(self):
        self.data = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00']

    def test_is_time_column_keyword(self):
        # Time column keyword without column number
        assert FITS_time.is_time_column_keyword('TRPOS') is False

        # Global time column keyword
        assert FITS_time.is_time_column_keyword('TIMESYS') is False

        # Valid time column keyword
        assert FITS_time.is_time_column_keyword('TRPOS12') is True

    def test_set_time(self):
        """
        Test that setting illegal time keywords fails.
        """
        ft = FITS_time()

        with pytest.raises(ValueError) as err:
            ft.set_global_time('TRPOS7', 8)
            assert 'Illegal global time keyword' in str(err.value)

        with pytest.raises(ValueError) as err:
            ft.set_column_override('TIMESYS', 8)
            assert 'Illegal column-specific time keyword' in str(err.value)

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_replace_time_table_loc(self, table_types):
        """
        Test all the unusual conditions for locations of Time columns in a Table.
        """
        t = table_types()
        t['a'] = Time(self.data, format='isot', scale='utc')
        t['b'] = Time(self.data, format='isot', scale='tt')

        # Check that vectorized location raises an exception
        t['a'].location = EarthLocation([1,2], [2,3], [3,4,])

        with pytest.raises(ValueError) as err:
            table, hdr = FITS_time.replace_time_table(t)
        assert 'Vectorized Location of Time Column' in str(err.value)

        # Check that multiple Time columns with different locations raise an exception
        t['a'].location = EarthLocation(1, 2, 3)
        t['b'].location = EarthLocation(2, 3, 4)

        with pytest.raises(ValueError) as err:
            table, hdr = FITS_time.replace_time_table(t)
            assert 'Multiple Time Columns with different geocentric' in str(err.value)

        # Check that Time column with no location specified will assume global location
        t['b'].location = None

        with catch_warnings() as w:
            table, hdr = FITS_time.replace_time_table(t)
            assert len(w) == 1
            assert str(w[0].message).startswith('Time Column "b" has no specified location,'
                                                ' but global Time Position is present')

        # Check that multiple Time columns with same location can be written
        t['b'].location = EarthLocation(1, 2, 3)

        with catch_warnings() as w:
            table, hdr = FITS_time.replace_time_table(t)
            assert len(w) == 0

        # Check compatibility of Time Scales and Reference Positions
        time_ref = FITS_time.TIME_SCALE_REF
        uncomp_scales = [scale for scale in time_ref if time_ref[scale] != 'TOPOCENTER']

        for scale in uncomp_scales:
            t.replace_column('a', getattr(t['a'], scale))
            with catch_warnings() as w:
                table, hdr = FITS_time.replace_time_table(t)
                assert len(w) == 1
                assert str(w[0].message).startswith('Earth Location "TOPOCENTER" for Time Column')

    @pytest.mark.parametrize('table_types', (Table, QTable))
    def test_replace_time_table_header(self, table_types):
        """
        Test the header returned by `replace_time_table` explicitly. 
        """
        t = table_types()
        t['a'] = Time(self.data, format='isot', scale='utc',
                      location=EarthLocation(-2446353.80003635,
                      4237209.07495215, 4077985.57220038, unit='m'))
        t['b'] = Time([1,2], format='cxcsec', scale='tt')

        ideal_col_hdr = {'OBSGEO-X' : t['a'].location.x.value,
                         'OBSGEO-Y' : t['a'].location.y.value,
                         'OBSGEO-Z' : t['a'].location.z.value,
                         'TCTYP1' : t['a'].scale.upper(),
                         'TRPOS1' : 'TOPOCENTER',
                         'TCTYP2' : t['b'].scale.upper()}
   
        table, hdr = FITS_time.replace_time_table(t)

        # Check global time keywords
        for key, value in FITS_time.GLOBAL_TIME_INFO.items():
             assert hdr[key] == value[0]
             assert hdr.comments[key] == value[1]
             hdr.remove(key)

        # Check column-specific(related) keywords
        for key, value in ideal_col_hdr.items():
            assert hdr[key] == value
            hdr.remove(key)

        assert len(hdr) == 2
