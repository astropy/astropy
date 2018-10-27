# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest
import numpy as np

from ..column import Column
from ..diff import (FITSDiff, HeaderDiff, ImageDataDiff, TableDataDiff,
                    HDUDiff)
from ..hdu import HDUList, PrimaryHDU, ImageHDU
from ..hdu.table import BinTableHDU
from ..header import Header

from ....tests.helper import catch_warnings
from ....utils.exceptions import AstropyDeprecationWarning
from ....io import fits

from . import FitsTestCase


class TestDiff(FitsTestCase):
    def test_identical_headers(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        assert HeaderDiff(ha, hb).identical
        assert HeaderDiff(ha.tostring(), hb.tostring()).identical

        with pytest.raises(TypeError):
            HeaderDiff(1, 2)

    def test_slightly_different_headers(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        assert not HeaderDiff(ha, hb).identical

    def test_common_keywords(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        hb['D'] = (5, 'Comment')
        assert HeaderDiff(ha, hb).common_keywords == ['A', 'B', 'C']

    def test_different_keyword_count(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        del hb['B']
        diff = HeaderDiff(ha, hb)
        assert not diff.identical
        assert diff.diff_keyword_count == (3, 2)

        # But make sure the common keywords are at least correct
        assert diff.common_keywords == ['A', 'C']

    def test_different_keywords(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        hb['D'] = (5, 'Comment')
        ha['E'] = (6, 'Comment')
        ha['F'] = (7, 'Comment')
        diff = HeaderDiff(ha, hb)
        assert not diff.identical
        assert diff.diff_keywords == (['E', 'F'], ['D'])

    def test_different_keyword_values(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        diff = HeaderDiff(ha, hb)
        assert not diff.identical
        assert diff.diff_keyword_values == {'C': [(3, 4)]}

    def test_different_keyword_comments(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3, 'comment 1')])
        hb = ha.copy()
        hb.comments['C'] = 'comment 2'
        diff = HeaderDiff(ha, hb)
        assert not diff.identical
        assert (diff.diff_keyword_comments ==
                {'C': [('comment 1', 'comment 2')]})

    def test_different_keyword_values_with_duplicate(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        ha.append(('C', 4))
        hb.append(('C', 5))
        diff = HeaderDiff(ha, hb)
        assert not diff.identical
        assert diff.diff_keyword_values == {'C': [None, (4, 5)]}

    def test_asymmetric_duplicate_keywords(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        ha.append(('A', 2, 'comment 1'))
        ha.append(('A', 3, 'comment 2'))
        hb.append(('B', 4, 'comment 3'))
        hb.append(('C', 5, 'comment 4'))
        diff = HeaderDiff(ha, hb)
        assert not diff.identical
        assert diff.diff_keyword_values == {}
        assert (diff.diff_duplicate_keywords ==
                {'A': (3, 1), 'B': (1, 2), 'C': (1, 2)})

        report = diff.report()
        assert ("Inconsistent duplicates of keyword 'A'     :\n"
                "  Occurs 3 time(s) in a, 1 times in (b)") in report

    def test_floating_point_rtol(self):
        ha = Header([('A', 1), ('B', 2.00001), ('C', 3.000001)])
        hb = ha.copy()
        hb['B'] = 2.00002
        hb['C'] = 3.000002
        diff = HeaderDiff(ha, hb)
        assert not diff.identical
        assert (diff.diff_keyword_values ==
                {'B': [(2.00001, 2.00002)], 'C': [(3.000001, 3.000002)]})
        diff = HeaderDiff(ha, hb, rtol=1e-6)
        assert not diff.identical
        assert diff.diff_keyword_values == {'B': [(2.00001, 2.00002)]}
        diff = HeaderDiff(ha, hb, rtol=1e-5)
        assert diff.identical

    def test_floating_point_atol(self):
        ha = Header([('A', 1), ('B', 1.0), ('C', 0.0)])
        hb = ha.copy()
        hb['B'] = 1.00001
        hb['C'] = 0.000001
        diff = HeaderDiff(ha, hb, rtol=1e-6)
        assert not diff.identical
        assert (diff.diff_keyword_values ==
                {'B': [(1.0, 1.00001)], 'C': [(0.0, 0.000001)]})
        diff = HeaderDiff(ha, hb, rtol=1e-5)
        assert not diff.identical
        assert (diff.diff_keyword_values ==
                {'C': [(0.0, 0.000001)]})
        diff = HeaderDiff(ha, hb, atol=1e-6)
        assert not diff.identical
        assert (diff.diff_keyword_values ==
                {'B': [(1.0, 1.00001)]})
        diff = HeaderDiff(ha, hb, atol=1e-5)  # strict inequality
        assert not diff.identical
        assert (diff.diff_keyword_values ==
                {'B': [(1.0, 1.00001)]})
        diff = HeaderDiff(ha, hb, rtol=1e-5, atol=1e-5)
        assert diff.identical
        diff = HeaderDiff(ha, hb, atol=1.1e-5)
        assert diff.identical
        diff = HeaderDiff(ha, hb, rtol=1e-6, atol=1e-6)
        assert not diff.identical

    def test_deprecation_tolerance(self):
        """Verify uses of tolerance and rtol.
        This test should be removed in the next astropy version."""

        ha = Header([('B', 1.0), ('C', 0.1)])
        hb = ha.copy()
        hb['B'] = 1.00001
        hb['C'] = 0.100001
        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            diff = HeaderDiff(ha, hb, tolerance=1e-6)
            assert warning_lines[0].category == AstropyDeprecationWarning
            assert (str(warning_lines[0].message) == '"tolerance" was '
                    'deprecated in version 2.0 and will be removed in a '
                    'future version. Use argument "rtol" instead.')
            assert (diff.diff_keyword_values == {'C': [(0.1, 0.100001)],
                                                 'B': [(1.0, 1.00001)]})
            assert not diff.identical

        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            # `rtol` is always ignored when `tolerance` is provided
            diff = HeaderDiff(ha, hb, rtol=1e-6, tolerance=1e-5)
            assert warning_lines[0].category == AstropyDeprecationWarning
            assert (str(warning_lines[0].message) == '"tolerance" was '
                    'deprecated in version 2.0 and will be removed in a '
                    'future version. Use argument "rtol" instead.')
            assert diff.identical

    def test_ignore_blanks(self):
        with fits.conf.set_temp('strip_header_whitespace', False):
            ha = Header([('A', 1), ('B', 2), ('C', 'A       ')])
            hb = ha.copy()
            hb['C'] = 'A'
            assert ha['C'] != hb['C']

            diff = HeaderDiff(ha, hb)
            # Trailing blanks are ignored by default
            assert diff.identical
            assert diff.diff_keyword_values == {}

            # Don't ignore blanks
            diff = HeaderDiff(ha, hb, ignore_blanks=False)
            assert not diff.identical
            assert diff.diff_keyword_values == {'C': [('A       ', 'A')]}

    def test_ignore_blank_cards(self):
        """Test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/152

        Ignore blank cards.
        """

        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = Header([('A', 1), ('', ''), ('B', 2), ('', ''), ('C', 3)])
        hc = ha.copy()
        hc.append()
        hc.append()

        # We now have a header with interleaved blanks, and a header with end
        # blanks, both of which should ignore the blanks
        assert HeaderDiff(ha, hb).identical
        assert HeaderDiff(ha, hc).identical
        assert HeaderDiff(hb, hc).identical

        assert not HeaderDiff(ha, hb, ignore_blank_cards=False).identical
        assert not HeaderDiff(ha, hc, ignore_blank_cards=False).identical

        # Both hb and hc have the same number of blank cards; since order is
        # currently ignored, these should still be identical even if blank
        # cards are not ignored
        assert HeaderDiff(hb, hc, ignore_blank_cards=False).identical

        hc.append()
        # But now there are different numbers of blanks, so they should not be
        # ignored:
        assert not HeaderDiff(hb, hc, ignore_blank_cards=False).identical

    def test_ignore_hdus(self):
        a = np.arange(100).reshape(10, 10)
        b = a.copy()
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        xa = np.array([(1.0, 1), (3.0, 4)], dtype=[('x', float), ('y', int)])
        xb = np.array([(1.0, 2), (3.0, 5)], dtype=[('x', float), ('y', int)])
        phdu = PrimaryHDU(header=ha)
        ihdua = ImageHDU(data=a, name='SCI')
        ihdub = ImageHDU(data=b, name='SCI')
        bhdu1 = BinTableHDU(data=xa, name='ASDF')
        bhdu2 = BinTableHDU(data=xb, name='ASDF')
        hdula = HDUList([phdu, ihdua, bhdu1])
        hdulb = HDUList([phdu, ihdub, bhdu2])

        # ASDF extension should be different
        diff = FITSDiff(hdula, hdulb)
        assert not diff.identical
        assert diff.diff_hdus[0][0] == 2

        # ASDF extension should be ignored
        diff = FITSDiff(hdula, hdulb, ignore_hdus=['ASDF'])
        assert diff.identical, diff.report()

        diff = FITSDiff(hdula, hdulb, ignore_hdus=['ASD*'])
        assert diff.identical, diff.report()

        # SCI extension should be different
        hdulb['SCI'].data += 1
        diff = FITSDiff(hdula, hdulb, ignore_hdus=['ASDF'])
        assert not diff.identical

        # SCI and ASDF extensions should be ignored
        diff = FITSDiff(hdula, hdulb, ignore_hdus=['SCI', 'ASDF'])
        assert diff.identical, diff.report()

        # All EXTVER of SCI should be ignored
        ihduc = ImageHDU(data=a, name='SCI', ver=2)
        hdulb.append(ihduc)
        diff = FITSDiff(hdula, hdulb, ignore_hdus=['SCI', 'ASDF'])
        assert not any(diff.diff_hdus), diff.report()
        assert any(diff.diff_hdu_count), diff.report()

    def test_ignore_keyword_values(self):
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['B'] = 4
        hb['C'] = 5
        diff = HeaderDiff(ha, hb, ignore_keywords=['*'])
        assert diff.identical
        diff = HeaderDiff(ha, hb, ignore_keywords=['B'])
        assert not diff.identical
        assert diff.diff_keyword_values == {'C': [(3, 5)]}

        report = diff.report()
        assert 'Keyword B        has different values' not in report
        assert 'Keyword C        has different values' in report

        # Test case-insensitivity
        diff = HeaderDiff(ha, hb, ignore_keywords=['b'])
        assert not diff.identical
        assert diff.diff_keyword_values == {'C': [(3, 5)]}

    def test_ignore_keyword_comments(self):
        ha = Header([('A', 1, 'A'), ('B', 2, 'B'), ('C', 3, 'C')])
        hb = ha.copy()
        hb.comments['B'] = 'D'
        hb.comments['C'] = 'E'
        diff = HeaderDiff(ha, hb, ignore_comments=['*'])
        assert diff.identical
        diff = HeaderDiff(ha, hb, ignore_comments=['B'])
        assert not diff.identical
        assert diff.diff_keyword_comments == {'C': [('C', 'E')]}

        report = diff.report()
        assert 'Keyword B        has different comments' not in report
        assert 'Keyword C        has different comments' in report

        # Test case-insensitivity
        diff = HeaderDiff(ha, hb, ignore_comments=['b'])
        assert not diff.identical
        assert diff.diff_keyword_comments == {'C': [('C', 'E')]}

    def test_trivial_identical_images(self):
        ia = np.arange(100).reshape(10, 10)
        ib = np.arange(100).reshape(10, 10)
        diff = ImageDataDiff(ia, ib)
        assert diff.identical
        assert diff.diff_total == 0

    def test_identical_within_relative_tolerance(self):
        ia = np.ones((10, 10)) - 0.00001
        ib = np.ones((10, 10)) - 0.00002
        diff = ImageDataDiff(ia, ib, rtol=1.0e-4)
        assert diff.identical
        assert diff.diff_total == 0

    def test_identical_within_absolute_tolerance(self):
        ia = np.zeros((10, 10)) - 0.00001
        ib = np.zeros((10, 10)) - 0.00002
        diff = ImageDataDiff(ia, ib, rtol=1.0e-4)
        assert not diff.identical
        assert diff.diff_total == 100
        diff = ImageDataDiff(ia, ib, atol=1.0e-4)
        assert diff.identical
        assert diff.diff_total == 0

    def test_identical_within_rtol_and_atol(self):
        ia = np.zeros((10, 10)) - 0.00001
        ib = np.zeros((10, 10)) - 0.00002
        diff = ImageDataDiff(ia, ib, rtol=1.0e-5, atol=1.0e-5)
        assert diff.identical
        assert diff.diff_total == 0

    def test_not_identical_within_rtol_and_atol(self):
        ia = np.zeros((10, 10)) - 0.00001
        ib = np.zeros((10, 10)) - 0.00002
        diff = ImageDataDiff(ia, ib, rtol=1.0e-5, atol=1.0e-6)
        assert not diff.identical
        assert diff.diff_total == 100

    def test_identical_comp_image_hdus(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/189

        For this test we mostly just care that comparing to compressed images
        does not crash, and returns the correct results.  Two compressed images
        will be considered identical if the decompressed data is the same.
        Obviously we test whether or not the same compression was used by
        looking for (or ignoring) header differences.
        """

        data = np.arange(100.0).reshape(10, 10)
        hdu = fits.CompImageHDU(data=data)
        hdu.writeto(self.temp('test.fits'))
        hdula = fits.open(self.temp('test.fits'))
        hdulb = fits.open(self.temp('test.fits'))
        diff = FITSDiff(hdula, hdulb)
        assert diff.identical

    def test_different_dimensions(self):
        ia = np.arange(100).reshape(10, 10)
        ib = np.arange(100) - 1

        # Although ib could be reshaped into the same dimensions, for now the
        # data is not compared anyways
        diff = ImageDataDiff(ia, ib)
        assert not diff.identical
        assert diff.diff_dimensions == ((10, 10), (100,))
        assert diff.diff_total == 0

        report = diff.report()
        assert 'Data dimensions differ' in report
        assert 'a: 10 x 10' in report
        assert 'b: 100' in report
        assert 'No further data comparison performed.'

    def test_different_pixels(self):
        ia = np.arange(100).reshape(10, 10)
        ib = np.arange(100).reshape(10, 10)
        ib[0, 0] = 10
        ib[5, 5] = 20
        diff = ImageDataDiff(ia, ib)
        assert not diff.identical
        assert diff.diff_dimensions == ()
        assert diff.diff_total == 2
        assert diff.diff_ratio == 0.02
        assert diff.diff_pixels == [((0, 0), (0, 10)), ((5, 5), (55, 20))]

    def test_identical_tables(self):
        c1 = Column('A', format='L', array=[True, False])
        c2 = Column('B', format='X', array=[[0], [1]])
        c3 = Column('C', format='4I', dim='(2, 2)',
                    array=[[0, 1, 2, 3], [4, 5, 6, 7]])
        c4 = Column('D', format='J', bscale=2.0, array=[0, 1])
        c5 = Column('E', format='A3', array=['abc', 'def'])
        c6 = Column('F', format='E', unit='m', array=[0.0, 1.0])
        c7 = Column('G', format='D', bzero=-0.1, array=[0.0, 1.0])
        c8 = Column('H', format='C', array=[0.0+1.0j, 2.0+3.0j])
        c9 = Column('I', format='M', array=[4.0+5.0j, 6.0+7.0j])
        c10 = Column('J', format='PI(2)', array=[[0, 1], [2, 3]])

        columns = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10]

        ta = BinTableHDU.from_columns(columns)
        tb = BinTableHDU.from_columns([c.copy() for c in columns])

        diff = TableDataDiff(ta.data, tb.data)
        assert diff.identical
        assert len(diff.common_columns) == 10
        assert diff.common_column_names == set('abcdefghij')
        assert diff.diff_ratio == 0
        assert diff.diff_total == 0

    def test_diff_empty_tables(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/178

        Ensure that diffing tables containing empty data doesn't crash.
        """

        c1 = Column('D', format='J')
        c2 = Column('E', format='J')
        thdu = BinTableHDU.from_columns([c1, c2], nrows=0)

        hdula = fits.HDUList([thdu])
        hdulb = fits.HDUList([thdu])

        diff = FITSDiff(hdula, hdulb)
        assert diff.identical

    def test_ignore_table_fields(self):
        c1 = Column('A', format='L', array=[True, False])
        c2 = Column('B', format='X', array=[[0], [1]])
        c3 = Column('C', format='4I', dim='(2, 2)',
                    array=[[0, 1, 2, 3], [4, 5, 6, 7]])

        c4 = Column('B', format='X', array=[[1], [0]])
        c5 = Column('C', format='4I', dim='(2, 2)',
                    array=[[1, 2, 3, 4], [5, 6, 7, 8]])

        ta = BinTableHDU.from_columns([c1, c2, c3])
        tb = BinTableHDU.from_columns([c1, c4, c5])

        diff = TableDataDiff(ta.data, tb.data, ignore_fields=['B', 'C'])
        assert diff.identical

        # The only common column should be c1
        assert len(diff.common_columns) == 1
        assert diff.common_column_names == {'a'}
        assert diff.diff_ratio == 0
        assert diff.diff_total == 0

    def test_different_table_field_names(self):
        ca = Column('A', format='L', array=[True, False])
        cb = Column('B', format='L', array=[True, False])
        cc = Column('C', format='L', array=[True, False])

        ta = BinTableHDU.from_columns([ca, cb])
        tb = BinTableHDU.from_columns([ca, cc])

        diff = TableDataDiff(ta.data, tb.data)

        assert not diff.identical
        assert len(diff.common_columns) == 1
        assert diff.common_column_names == {'a'}
        assert diff.diff_column_names == (['B'], ['C'])
        assert diff.diff_ratio == 0
        assert diff.diff_total == 0

        report = diff.report()
        assert 'Extra column B of format L in a' in report
        assert 'Extra column C of format L in b' in report

    def test_different_table_field_counts(self):
        """
        Test tables with some common columns, but different number of columns
        overall.
        """

        ca = Column('A', format='L', array=[True, False])
        cb = Column('B', format='L', array=[True, False])
        cc = Column('C', format='L', array=[True, False])

        ta = BinTableHDU.from_columns([cb])
        tb = BinTableHDU.from_columns([ca, cb, cc])

        diff = TableDataDiff(ta.data, tb.data)

        assert not diff.identical
        assert diff.diff_column_count == (1, 3)
        assert len(diff.common_columns) == 1
        assert diff.common_column_names == {'b'}
        assert diff.diff_column_names == ([], ['A', 'C'])
        assert diff.diff_ratio == 0
        assert diff.diff_total == 0

        report = diff.report()
        assert ' Tables have different number of columns:' in report
        assert '  a: 1\n  b: 3' in report

    def test_different_table_rows(self):
        """
        Test tables that are otherwise identical but one has more rows than the
        other.
        """

        ca1 = Column('A', format='L', array=[True, False])
        cb1 = Column('B', format='L', array=[True, False])
        ca2 = Column('A', format='L', array=[True, False, True])
        cb2 = Column('B', format='L', array=[True, False, True])

        ta = BinTableHDU.from_columns([ca1, cb1])
        tb = BinTableHDU.from_columns([ca2, cb2])

        diff = TableDataDiff(ta.data, tb.data)

        assert not diff.identical
        assert diff.diff_column_count == ()
        assert len(diff.common_columns) == 2
        assert diff.diff_rows == (2, 3)
        assert diff.diff_values == []

        report = diff.report()

        assert 'Table rows differ' in report
        assert 'a: 2' in report
        assert 'b: 3' in report
        assert 'No further data comparison performed.'

    def test_different_table_data(self):
        """
        Test diffing table data on columns of several different data formats
        and dimensions.
        """

        ca1 = Column('A', format='L', array=[True, False])
        ca2 = Column('B', format='X', array=[[0], [1]])
        ca3 = Column('C', format='4I', dim='(2, 2)',
                     array=[[0, 1, 2, 3], [4, 5, 6, 7]])
        ca4 = Column('D', format='J', bscale=2.0, array=[0.0, 2.0])
        ca5 = Column('E', format='A3', array=['abc', 'def'])
        ca6 = Column('F', format='E', unit='m', array=[0.0, 1.0])
        ca7 = Column('G', format='D', bzero=-0.1, array=[0.0, 1.0])
        ca8 = Column('H', format='C', array=[0.0+1.0j, 2.0+3.0j])
        ca9 = Column('I', format='M', array=[4.0+5.0j, 6.0+7.0j])
        ca10 = Column('J', format='PI(2)', array=[[0, 1], [2, 3]])

        cb1 = Column('A', format='L', array=[False, False])
        cb2 = Column('B', format='X', array=[[0], [0]])
        cb3 = Column('C', format='4I', dim='(2, 2)',
                     array=[[0, 1, 2, 3], [5, 6, 7, 8]])
        cb4 = Column('D', format='J', bscale=2.0, array=[2.0, 2.0])
        cb5 = Column('E', format='A3', array=['abc', 'ghi'])
        cb6 = Column('F', format='E', unit='m', array=[1.0, 2.0])
        cb7 = Column('G', format='D', bzero=-0.1, array=[2.0, 3.0])
        cb8 = Column('H', format='C', array=[1.0+1.0j, 2.0+3.0j])
        cb9 = Column('I', format='M', array=[5.0+5.0j, 6.0+7.0j])
        cb10 = Column('J', format='PI(2)', array=[[1, 2], [3, 4]])

        ta = BinTableHDU.from_columns([ca1, ca2, ca3, ca4, ca5, ca6, ca7,
                                       ca8, ca9, ca10])
        tb = BinTableHDU.from_columns([cb1, cb2, cb3, cb4, cb5, cb6, cb7,
                                       cb8, cb9, cb10])

        diff = TableDataDiff(ta.data, tb.data, numdiffs=20)
        assert not diff.identical
        # The column definitions are the same, but not the column values
        assert diff.diff_columns == ()
        assert diff.diff_values[0] == (('A', 0), (True, False))
        assert diff.diff_values[1] == (('B', 1), ([1], [0]))
        assert diff.diff_values[2][0] == ('C', 1)
        assert (diff.diff_values[2][1][0] == [[4, 5], [6, 7]]).all()
        assert (diff.diff_values[2][1][1] == [[5, 6], [7, 8]]).all()
        assert diff.diff_values[3] == (('D', 0), (0, 2.0))
        assert diff.diff_values[4] == (('E', 1), ('def', 'ghi'))
        assert diff.diff_values[5] == (('F', 0), (0.0, 1.0))
        assert diff.diff_values[6] == (('F', 1), (1.0, 2.0))
        assert diff.diff_values[7] == (('G', 0), (0.0, 2.0))
        assert diff.diff_values[8] == (('G', 1), (1.0, 3.0))
        assert diff.diff_values[9] == (('H', 0), (0.0+1.0j, 1.0+1.0j))
        assert diff.diff_values[10] == (('I', 0), (4.0+5.0j, 5.0+5.0j))
        assert diff.diff_values[11][0] == ('J', 0)
        assert (diff.diff_values[11][1][0] == [0, 1]).all()
        assert (diff.diff_values[11][1][1] == [1, 2]).all()
        assert diff.diff_values[12][0] == ('J', 1)
        assert (diff.diff_values[12][1][0] == [2, 3]).all()
        assert (diff.diff_values[12][1][1] == [3, 4]).all()

        assert diff.diff_total == 13
        assert diff.diff_ratio == 0.65

        report = diff.report()
        assert ('Column A data differs in row 0:\n'
                '    a> True\n'
                '    b> False') in report
        assert ('...and at 1 more indices.\n'
                ' Column D data differs in row 0:') in report
        assert ('13 different table data element(s) found (65.00% different)'
                in report)
        assert report.count('more indices') == 1

    def test_identical_files_basic(self):
        """Test identicality of two simple, extensionless files."""

        a = np.arange(100).reshape(10, 10)
        hdu = PrimaryHDU(data=a)
        hdu.writeto(self.temp('testa.fits'))
        hdu.writeto(self.temp('testb.fits'))
        diff = FITSDiff(self.temp('testa.fits'), self.temp('testb.fits'))
        assert diff.identical

        report = diff.report()
        # Primary HDUs should contain no differences
        assert 'Primary HDU' not in report
        assert 'Extension HDU' not in report
        assert 'No differences found.' in report

        a = np.arange(10)
        ehdu = ImageHDU(data=a)
        diff = HDUDiff(ehdu, ehdu)
        assert diff.identical
        report = diff.report()
        assert 'No differences found.' in report

    def test_partially_identical_files1(self):
        """
        Test files that have some identical HDUs but a different extension
        count.
        """

        a = np.arange(100).reshape(10, 10)
        phdu = PrimaryHDU(data=a)
        ehdu = ImageHDU(data=a)
        hdula = HDUList([phdu, ehdu])
        hdulb = HDUList([phdu, ehdu, ehdu])
        diff = FITSDiff(hdula, hdulb)
        assert not diff.identical
        assert diff.diff_hdu_count == (2, 3)

        # diff_hdus should be empty, since the third extension in hdulb
        # has nothing to compare against
        assert diff.diff_hdus == []

        report = diff.report()
        assert 'Files contain different numbers of HDUs' in report
        assert 'a: 2\n b: 3' in report
        assert 'No differences found between common HDUs' in report

    def test_partially_identical_files2(self):
        """
        Test files that have some identical HDUs but one different HDU.
        """

        a = np.arange(100).reshape(10, 10)
        phdu = PrimaryHDU(data=a)
        ehdu = ImageHDU(data=a)
        ehdu2 = ImageHDU(data=(a + 1))
        hdula = HDUList([phdu, ehdu, ehdu])
        hdulb = HDUList([phdu, ehdu2, ehdu])
        diff = FITSDiff(hdula, hdulb)

        assert not diff.identical
        assert diff.diff_hdu_count == ()
        assert len(diff.diff_hdus) == 1
        assert diff.diff_hdus[0][0] == 1

        hdudiff = diff.diff_hdus[0][1]
        assert not hdudiff.identical
        assert hdudiff.diff_extnames == ()
        assert hdudiff.diff_extvers == ()
        assert hdudiff.diff_extension_types == ()
        assert hdudiff.diff_headers.identical
        assert hdudiff.diff_data is not None

        datadiff = hdudiff.diff_data
        assert isinstance(datadiff, ImageDataDiff)
        assert not datadiff.identical
        assert datadiff.diff_dimensions == ()
        assert (datadiff.diff_pixels ==
                [((0, y), (y, y + 1)) for y in range(10)])
        assert datadiff.diff_ratio == 1.0
        assert datadiff.diff_total == 100

        report = diff.report()
        # Primary HDU and 2nd extension HDU should have no differences
        assert 'Primary HDU' not in report
        assert 'Extension HDU 2' not in report
        assert 'Extension HDU 1' in report

        assert 'Headers contain differences' not in report
        assert 'Data contains differences' in report
        for y in range(10):
            assert 'Data differs at [{}, 1]'.format(y + 1) in report
        assert '100 different pixels found (100.00% different).' in report

    def test_partially_identical_files3(self):
        """
        Test files that have some identical HDUs but a different extension
        name.
        """

        phdu = PrimaryHDU()
        ehdu = ImageHDU(name='FOO')
        hdula = HDUList([phdu, ehdu])
        ehdu = BinTableHDU(name='BAR')
        ehdu.header['EXTVER'] = 2
        ehdu.header['EXTLEVEL'] = 3
        hdulb = HDUList([phdu, ehdu])
        diff = FITSDiff(hdula, hdulb)
        assert not diff.identical

        assert diff.diff_hdus[0][0] == 1

        hdu_diff = diff.diff_hdus[0][1]
        assert hdu_diff.diff_extension_types == ('IMAGE', 'BINTABLE')
        assert hdu_diff.diff_extnames == ('FOO', 'BAR')
        assert hdu_diff.diff_extvers == (1, 2)
        assert hdu_diff.diff_extlevels == (1, 3)

        report = diff.report()
        assert 'Extension types differ' in report
        assert 'a: IMAGE\n    b: BINTABLE' in report
        assert 'Extension names differ' in report
        assert 'a: FOO\n    b: BAR' in report
        assert 'Extension versions differ' in report
        assert 'a: 1\n    b: 2' in report
        assert 'Extension levels differ' in report
        assert 'a: 1\n    b: 2' in report

    def test_diff_nans(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/204
        """

        # First test some arrays that should be equivalent....
        arr = np.empty((10, 10), dtype=np.float64)
        arr[:5] = 1.0
        arr[5:] = np.nan
        arr2 = arr.copy()

        table = np.rec.array([(1.0, 2.0), (3.0, np.nan), (np.nan, np.nan)],
                             names=['cola', 'colb']).view(fits.FITS_rec)
        table2 = table.copy()

        assert ImageDataDiff(arr, arr2).identical
        assert TableDataDiff(table, table2).identical

        # Now let's introduce some differences, where there are nans and where
        # there are not nans
        arr2[0][0] = 2.0
        arr2[5][0] = 2.0
        table2[0][0] = 2.0
        table2[1][1] = 2.0

        diff = ImageDataDiff(arr, arr2)
        assert not diff.identical
        assert diff.diff_pixels[0] == ((0, 0), (1.0, 2.0))
        assert diff.diff_pixels[1][0] == (5, 0)
        assert np.isnan(diff.diff_pixels[1][1][0])
        assert diff.diff_pixels[1][1][1] == 2.0

        diff = TableDataDiff(table, table2)
        assert not diff.identical
        assert diff.diff_values[0] == (('cola', 0), (1.0, 2.0))
        assert diff.diff_values[1][0] == ('colb', 1)
        assert np.isnan(diff.diff_values[1][1][0])
        assert diff.diff_values[1][1][1] == 2.0

    def test_file_output_from_path_string(self):
        outpath = self.temp('diff_output.txt')
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        diffobj = HeaderDiff(ha, hb)
        diffobj.report(fileobj=outpath)
        report_as_string = diffobj.report()
        with open(outpath) as fout:
            assert fout.read() == report_as_string

    def test_file_output_overwrite_safety(self):
        outpath = self.temp('diff_output.txt')
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        diffobj = HeaderDiff(ha, hb)
        diffobj.report(fileobj=outpath)

        with pytest.raises(OSError):
            diffobj.report(fileobj=outpath)

    def test_file_output_overwrite_success(self):
        outpath = self.temp('diff_output.txt')
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        diffobj = HeaderDiff(ha, hb)
        diffobj.report(fileobj=outpath)
        report_as_string = diffobj.report()
        diffobj.report(fileobj=outpath, overwrite=True)
        with open(outpath) as fout:
            assert fout.read() == report_as_string, (
                "overwritten output file is not identical to report string")

    def test_file_output_overwrite_vs_clobber(self):
        """Verify uses of clobber and overwrite."""

        outpath = self.temp('diff_output.txt')
        ha = Header([('A', 1), ('B', 2), ('C', 3)])
        hb = ha.copy()
        hb['C'] = 4
        diffobj = HeaderDiff(ha, hb)
        diffobj.report(fileobj=outpath)
        with catch_warnings(AstropyDeprecationWarning) as warning_lines:
            diffobj.report(fileobj=outpath, clobber=True)
            assert warning_lines[0].category == AstropyDeprecationWarning
            assert (str(warning_lines[0].message) == '"clobber" was '
                    'deprecated in version 2.0 and will be removed in a '
                    'future version. Use argument "overwrite" instead.')
