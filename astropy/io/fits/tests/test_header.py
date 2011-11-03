import pyfits

from pyfits.tests import PyfitsTestCase

from nose.tools import assert_equal, assert_false, assert_raises, assert_true


class TestRecordValuedKeywordCards(PyfitsTestCase):
    """Tests for handling of record-valued keyword cards as used by the FITS
    WCS Paper IV proposal.

    These tests are derived primarily from the release notes for PyFITS 1.4 (in
    which this feature was first introduced.
    """

    def setup(self):
        super(TestRecordValuedKeywordCards, self).setup()
        self._test_header = pyfits.Header()
        self._test_header.update('DP1', 'NAXIS: 2')
        self._test_header.update('DP1', 'AXIS.1: 1')
        self._test_header.update('DP1', 'AXIS.2: 2')
        self._test_header.update('DP1', 'NAUX: 2')
        self._test_header.update('DP1', 'AUX.1.COEFF.0: 0')
        self._test_header.update('DP1', 'AUX.1.POWER.0: 1')
        self._test_header.update('DP1', 'AUX.1.COEFF.1: 0.00048828125')
        self._test_header.update('DP1', 'AUX.1.POWER.1: 1')

    def test_field_specifier_case_senstivity(self):
        """The keyword portion of an RVKC should still be case-insensitive, but
        the the field-specifier portion should be case-sensitive.
        """

        header = pyfits.Header()
        header.update('abc.def', 1)
        header.update('abc.DEF', 2)
        assert_equal(header['abc.def'], 1)
        assert_equal(header['ABC.def'], 1)
        assert_equal(header['aBc.def'], 1)
        assert_equal(header['ABC.DEF'], 2)
        assert_false('ABC.dEf' in header)

    def test_get_rvkc_by_index(self):
        """Returning a RVKC from a header via index lookup should return the
        entire string value of the card, including the field-specifier.
        """

        assert_equal(self._test_header[0], 'NAXIS: 2')
        assert_equal(self._test_header[1], 'AXIS.1: 1')

    def test_get_rvkc_by_keyword(self):
        """Returning a RVKC just via the keyword name should return the entire
        string value of the first card with that keyword.
        """

        assert_equal(self._test_header['DP1'], 'NAXIS: 2')

    def test_get_rvkc_by_keyword_and_field_specifier(self):
        """Returning a RVKC via the full keyword/field-specifier combination
        should return the floating point value associated with the RVKC.
        """

        assert_equal(self._test_header['DP1.NAXIS'], 2.0)
        assert_true(isinstance(self._test_header['DP1.NAXIS'], float))
        assert_equal(self._test_header['DP1.AUX.1.COEFF.1'], 0.00048828125)

    def test_access_nonexistent_rvkc(self):
        """Accessing a nonexistent RVKC should raise an IndexError for
        index-based lookup, or a KeyError for keyword lookup (like a normal
        card).
        """

        assert_raises(IndexError, lambda x: self._test_header[x], 8)
        assert_raises(KeyError, lambda k: self._test_header[k], 'DP1.AXIS.3')

    def test_update_rvkc(self):
        """A RVKC can be updated either via index or keyword access."""

        self._test_header[0] = 3
        assert_equal(self._test_header['DP1.NAXIS'], 3.0)
        assert_true(isinstance(self._test_header['DP1.NAXIS'], float))
        assert_equal(self._test_header[0], 'NAXIS: 3')

        self._test_header['DP1.AXIS.1'] = 1.1
        assert_equal(self._test_header['DP1.AXIS.1'], 1.1)
        assert_equal(self._test_header[1], 'AXIS.1: 1.1')

    def test_rvkc_insert_after(self):
        """It should be possible to insert a new RVKC after an existing one
        specified by the full keyword/field-specifier combination."""

        self._test_header.update('DP1', 'AXIS.3: 1', 'a comment',
                                 after='DP1.AXIS.2')
        assert_equal(self._test_header[3], 'AXIS.3: 1')
        assert_equal(self._test_header['DP1.AXIS.3'], 1)

    def test_rvkc_delete(self):
        """Deleting a RVKC should work as with a normal card by using the full
        keyword/field-spcifier combination.
        """

        del self._test_header['DP1.AXIS.1']
        assert_equal(len(self._test_header), 7)
        assert_equal(self._test_header[0], 'NAXIS: 2')
        assert_equal(self._test_header[1], 'AXIS.2: 2')

    def test_pattern_matching_keys(self):
        """Test the keyword filter strings with RVKCs."""

        cl = self._test_header['DP1.AXIS.*']
        assert_true(isinstance(cl, pyfits.CardList))
        assert_equal(
            [l.strip() for l in str(cl).splitlines()],
            ["DP1     = 'AXIS.1: 1'",
             "DP1     = 'AXIS.2: 2'"])

        cl = self._test_header['DP1.N*']
        assert_equal(
            [l.strip() for l in str(cl).splitlines()],
            ["DP1     = 'NAXIS: 2'",
             "DP1     = 'NAUX: 2'"])

        cl = self._test_header['DP1.AUX...']
        assert_equal(
            [l.strip() for l in str(cl).splitlines()],
            ["DP1     = 'AUX.1.COEFF.0: 0'",
             "DP1     = 'AUX.1.POWER.0: 1'",
             "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
             "DP1     = 'AUX.1.POWER.1: 1'"])

        cl = self._test_header['DP?.NAXIS']
        assert_equal(
            [l.strip() for l in str(cl).splitlines()],
            ["DP1     = 'NAXIS: 2'"])

        cl = self._test_header['DP1.A*S.*']
        assert_equal(
            [l.strip() for l in str(cl).splitlines()],
            ["DP1     = 'AXIS.1: 1'",
             "DP1     = 'AXIS.2: 2'"])

    def test_pattern_matching_key_deletion(self):
        """Deletion by filter strings should work."""

        del self._test_header['DP1.A*...']
        assert_equal(len(self._test_header), 2)
        assert_equal(self._test_header[0], 'NAXIS: 2')
        assert_equal(self._test_header[1], 'NAUX: 2')

    def test_successive_pattern_matching(self):
        """A card list returned via a filter string should be further
        filterable."""

        cl = self._test_header['DP1.A*...']
        assert_equal(
            [l.strip() for l in str(cl).splitlines()],
            ["DP1     = 'AXIS.1: 1'",
             "DP1     = 'AXIS.2: 2'",
             "DP1     = 'AUX.1.COEFF.0: 0'",
             "DP1     = 'AUX.1.POWER.0: 1'",
             "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
             "DP1     = 'AUX.1.POWER.1: 1'"])

        cl2 = cl['*.*AUX...']
        assert_equal(
            [l.strip() for l in str(cl2).splitlines()],
            ["DP1     = 'AUX.1.COEFF.0: 0'",
             "DP1     = 'AUX.1.POWER.0: 1'",
             "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
             "DP1     = 'AUX.1.POWER.1: 1'"])

    def test_rvkc_in_cardlist_keys(self):
        """The CardList.keys() method should return full keyword/field-spec
        values for RVKCs.
        """

        cl = self._test_header['DP1.AXIS.*']
        assert_equal(cl.keys(), ['DP1.AXIS.1', 'DP1.AXIS.2'])

    def test_rvkc_in_cardlist_values(self):
        """The CardList.values() method should return the values of all RVKCs
        as floating point values.
        """

        cl = self._test_header['DP1.AXIS.*']
        assert_equal(cl.values(), [1.0, 2.0])

    def test_rvkc_cardlist_indexing(self):
        """RVKC should be retrievable from CardLists using standard index or
        keyword-based lookup.
        """

        cl = self._test_header['DP1.AXIS.*']
        assert_equal(str(cl[0]).strip(), "DP1     = 'AXIS.1: 1'")
        assert_equal(str(cl['DP1.AXIS.2']).strip(), "DP1     = 'AXIS.2: 2'")

    def test_rvkc_value_attribute(self):
        """Individual card values should be accessible by the .value attribute
        (which should return a float).
        """

        cl = self._test_header['DP1.AXIS.*']
        assert_equal(cl[0].value, 1.0)
        assert_true(isinstance(cl[0].value, float))

    def test_rvkc_cardlist_deletion(self):
        """Modifying RVKCs in a CardList should be reflected in any Header that
        those cards belong to, but deleting a RVKC from a CardList should not
        remove that card from the Header.
        """

        cl = self._test_header['DP1.AXIS.*']
        cl[0].value = 4.0
        assert_equal(self._test_header['DP1.AXIS.1'], 4.0)
        del cl[0]
        assert_raises(KeyError, lambda k: cl[k], 'DP1.AXIS.1')
        assert_equal(self._test_header['DP1.AXIS.1'], 4.0)

    def test_rvkc_constructor(self):
        """Test direct creation of RVKC objects."""

        c1 = pyfits.RecordValuedKeywordCard('DP1', 'NAXIS: 2',
                                            'Number of variables')
        c2 = pyfits.RecordValuedKeywordCard('DP1.AXIS.1', 1.0, 'Axis number')

        assert_equal(c1.key, 'DP1.NAXIS')
        assert_equal(c1.value, 2.0)
        assert_equal(c1.comment, 'Number of variables')
        assert_equal(c1.field_specifier, 'NAXIS')
        assert_equal(c2.key, 'DP1.AXIS.1')
        assert_equal(c2.value, 1.0)
        assert_equal(c2.comment, 'Axis number')
        assert_equal(c2.field_specifier, 'AXIS.1')

        # RVKCs created with this constructor should verify without any
        # problems
        c1.verify('exception')
        c2.verify('exception')

    def test_rvkc_fromstring(self):
        """Test creation of RVKC from their string representation."""

        c1 = pyfits.RecordValuedKeywordCard().fromstring(
            "DP1     = 'NAXIS: 2' / Number of independent variables")
        c2 = pyfits.RecordValuedKeywordCard().fromstring(
            "DP1     = 'AXIS.1: X' / Axis number")

        assert_equal(str(c1).strip(),
                     "DP1     = 'NAXIS: 2' / Number of independent variables")
        assert_equal(str(c2).strip(),
                     "DP1     = 'AXIS.1: X' / Axis number")
        # Since c2's value is wrong for a RVKC it should be a normal card
        assert_false(isinstance(c2, pyfits.RecordValuedKeywordCard))
