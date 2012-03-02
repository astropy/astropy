# Licensed under a 3-clause BSD style license - see PYFITS.rst

import itertools
import warnings

from ....io import fits
from ....tests.helper import pytest

from . import FitsTestCase
from .util import ignore_warnings
from ..card import _pad


class TestOldApiHeaderFunctions(FitsTestCase):
    """
    Tests that specifically use attributes and methods from the old
    Header/CardList API from PyFITS 3.0 and prior.

    This tests backward compatibility support for those interfaces.
    """

    def test_ascardimage_verifies_the_comment_string_to_be_ascii_text(self):
        # the ascardimage() verifies the comment string to be ASCII text
        c = fits.Card.fromstring('abc     = +  2.1   e + 12 / abcde\0')
        pytest.raises(Exception, c.ascardimage)

    def test_rename_key(self):
        """Test backwards compatibility support for Header.rename_key()"""
        header = fits.Header([('A', 'B', 'C'), ('D', 'E', 'F')])
        header.rename_key('A', 'B')
        assert 'A' not in header
        assert 'B' in header
        assert header[0] == 'B'
        assert header['B'] == 'B'
        assert header.comments['B'] == 'C'

    def test_add_commentary(self):
        header = fits.Header([('A', 'B', 'C'), ('HISTORY', 1),
                              ('HISTORY', 2), ('HISTORY', 3), ('', '', ''),
                              ('', '', '')])
        header.add_history(4)
        # One of the blanks should get used, so the length shouldn't change
        assert len(header) == 6
        assert header.cards[4].value == 4
        assert header['HISTORY'] == [1, 2, 3, 4]

        header.add_history(0, after='A')
        assert len(header) == 6
        assert header.cards[1].value == 0
        assert header['HISTORY'] == [0, 1, 2, 3, 4]

        header = fits.Header([('A', 'B', 'C'), ('', 1), ('', 2), ('', 3),
                              ('', '', ''), ('', '', '')])
        header.add_blank(4)
        # This time a new blank should be added, and the existing blanks don't
        # get used... (though this is really kinda sketchy--there's a
        # distinction between truly blank cards, and cards with blank keywords
        # that isn't currently made int he code)
        assert len(header) == 7
        assert header.cards[6].value == 4
        assert header[''] == [1, 2, 3, '', '', 4]

        header.add_blank(0, after='A')
        assert len(header) == 8
        assert header.cards[1].value == 0
        assert header[''] == [0, 1, 2, 3, '', '', 4]

    def test_has_key(self):
        header = fits.Header([('A', 'B', 'C'), ('D', 'E', 'F')])
        assert header.has_key('A')
        assert header.has_key('D')
        assert not header.has_key('C')

    def test_totxtfile(self):
        hdul = fits.open(self.data('test0.fits'))
        hdul[0].header.toTxtFile(self.temp('header.txt'))
        hdu = fits.ImageHDU()
        hdu.header.update('MYKEY', 'FOO', 'BAR')
        hdu.header.fromTxtFile(self.temp('header.txt'), replace=True)
        assert len(hdul[0].header.ascard) == len(hdu.header.ascard)
        assert not hdu.header.has_key('MYKEY')
        assert not hdu.header.has_key('EXTENSION')
        assert hdu.header.has_key('SIMPLE')

        # Write the hdu out and read it back in again--it should be recognized
        # as a PrimaryHDU
        hdu.writeto(self.temp('test.fits'), output_verify='ignore')
        assert (isinstance(fits.open(self.temp('test.fits'))[0],
                           fits.PrimaryHDU))

        hdu = fits.ImageHDU()
        hdu.header.update('MYKEY', 'FOO', 'BAR')
        hdu.header.fromTxtFile(self.temp('header.txt'))
        # hdu.header should have MYKEY keyword, and also adds PCOUNT and
        # GCOUNT, giving it 3 more keywords in total than the original
        assert len(hdul[0].header.ascard) == len(hdu.header.ascard) - 3
        assert hdu.header.has_key('MYKEY')
        assert not hdu.header.has_key('EXTENSION')
        assert hdu.header.has_key('SIMPLE')

        with ignore_warnings():
            hdu.writeto(self.temp('test.fits'), output_verify='ignore',
                        clobber=True)
        hdul2 = fits.open(self.temp('test.fits'))
        assert len(hdul2) == 2
        assert hdul2[1].header.has_key('MYKEY')

    def test_update_comment(self):
        hdul = fits.open(self.data('arange.fits'))
        hdul[0].header.update('FOO', 'BAR', 'BAZ')
        assert hdul[0].header['FOO'] == 'BAR'
        assert hdul[0].header.ascard['FOO'].comment == 'BAZ'

        hdul.writeto(self.temp('test.fits'))

        hdul = fits.open(self.temp('test.fits'), mode='update')
        hdul[0].header.ascard['FOO'].comment = 'QUX'
        hdul.close()

        hdul = fits.open(self.temp('test.fits'))
        assert hdul[0].header.ascard['FOO'].comment == 'QUX'

    def test_long_commentary_card(self):
        # Another version of this test using new API methods is found in
        # TestHeaderFunctions
        header = fits.Header()
        header.update('FOO', 'BAR')
        header.update('BAZ', 'QUX')
        longval = 'ABC' * 30
        header.add_history(longval)
        header.update('FRED', 'BARNEY')
        header.add_history(longval)

        assert len(header.ascard) == 7
        assert header.ascard[2].key == 'FRED'
        assert str(header.cards[3]) == 'HISTORY ' + longval[:72]
        assert str(header.cards[4]).rstrip() == 'HISTORY ' + longval[72:]

        header.add_history(longval, after='FOO')
        assert len(header.ascard) == 9
        assert str(header.cards[1]) == 'HISTORY ' + longval[:72]
        assert str(header.cards[2]).rstrip() == 'HISTORY ' + longval[72:]

    def test_wildcard_slice(self):
        """Test selecting a subsection of a header via wildcard matching."""

        header = fits.Header()
        header.update('ABC', 0)
        header.update('DEF', 1)
        header.update('ABD', 2)
        cards = header.ascard['AB*']
        assert len(cards) == 2
        assert cards[0].value == 0
        assert cards[1].value == 2


class TestHeaderFunctions(FitsTestCase):
    """Test Header and Card objects."""

    def test_card_constructor_default_args(self):
        """Test Card constructor with default argument values."""

        c = fits.Card()
        assert '' == c.key

    def test_string_value_card(self):
        """Test Card constructor with string value"""

        c = fits.Card('abc', '<8 ch')
        assert str(c) == _pad("ABC     = '<8 ch   '")
        c = fits.Card('nullstr', '')
        assert str(c) == _pad("NULLSTR = ''")

    def test_boolean_value_card(self):
        """Test Card constructor with boolean value"""

        c = fits.Card("abc", True)
        assert str(c) == _pad("ABC     =                    T")

        c = fits.Card.fromstring('abc     = F')
        assert c.value == False

    def test_long_integer_value_card(self):
        """Test Card constructor with long integer value"""

        c = fits.Card('long_int', -467374636747637647347374734737437)
        assert str(c) == _pad("LONG_INT= -467374636747637647347374734737437")

    def test_floating_point_value_card(self):
        """Test Card constructor with floating point value"""

        c = fits.Card('floatnum', -467374636747637647347374734737437.)

        if (str(c) != _pad("FLOATNUM= -4.6737463674763E+32") and
            str(c) != _pad("FLOATNUM= -4.6737463674763E+032")):
            assert str(c) == _pad("FLOATNUM= -4.6737463674763E+32")

    def test_complex_value_card(self):
        """Test Card constructor with complex value"""

        c = fits.Card('abc',
                      (1.2345377437887837487e88 + 6324767364763746367e-33j))
        f1 = _pad("ABC     = (1.23453774378878E+88, 6.32476736476374E-15)")
        f2 = _pad("ABC     = (1.2345377437887E+088, 6.3247673647637E-015)")
        f3 = _pad("ABC     = (1.23453774378878E+88, 6.32476736476374E-15)")
        if str(c) != f1 and str(c) != f2:
            assert str(c) == f3

    def test_card_image_constructed_too_long(self):
        """Test that over-long cards truncate the comment"""

        # card image constructed from key/value/comment is too long
        # (non-string value)
        with ignore_warnings():
            c = fits.Card('abc', 9, 'abcde' * 20)
            assert (str(c) ==
                    "ABC     =                    9 "
                    "/ abcdeabcdeabcdeabcdeabcdeabcdeabcdeabcdeabcdeab")
            c = fits.Card('abc', 'a' * 68, 'abcdefg')
            assert str(c) == "ABC     = '%s'" % ('a' * 68)

    def test_constructor_filter_illegal_data_structures(self):
        """Test that Card constructor raises exceptions on bad arguments"""

        pytest.raises(ValueError, fits.Card, ('abc',), {'value': (2, 3)})
        pytest.raises(ValueError, fits.Card, 'key', [], 'comment')

    def test_keyword_too_long(self):
        """Test that long Card keywords are allowed, but with a warning"""

        with warnings.catch_warnings():
            warnings.simplefilter('error')
            pytest.raises(UserWarning, fits.Card, 'abcdefghi', 'long')

    def test_illegal_characters_in_key(self):
        """
        Test that Card constructor disallows illegal characters in the keyword
        """

        pytest.raises(ValueError, fits.Card, 'abc+', 9)

    def test_commentary_cards(self):
        # commentary cards
        val = "A commentary card's value has no quotes around it."
        c = fits.Card("history", val)
        assert str(c) == _pad('HISTORY ' + val)
        val = "A commentary card has no comment."
        c = fits.Card("comment", val, "comment")
        assert str(c) == _pad('COMMENT ' + val)

    def test_commentary_card_created_by_fromstring(self):
        # commentary card created by fromstring()
        c = fits.Card.fromstring(
            "COMMENT card has no comments. "
            "/ text after slash is still part of the value.")
        assert (c.value == 'card has no comments. '
                           '/ text after slash is still part of the value.')
        assert c.comment == ''

    def test_commentary_card_will_not_parse_numerical_value(self):
        # commentary card will not parse the numerical value
        c = fits.Card.fromstring("history  (1, 2)")
        assert str(c) == _pad("HISTORY  (1, 2)")

    def test_equal_sign_after_column8(self):
        # equal sign after column 8 of a commentary card will be part ofthe
        # string value
        c = fits.Card.fromstring("history =   (1, 2)")
        assert str(c) == _pad("HISTORY =   (1, 2)")

    def test_blank_keyword(self):
        c = fits.Card('', '       / EXPOSURE INFORMATION')
        assert str(c) == _pad('               / EXPOSURE INFORMATION')
        c = fits.Card.fromstring(str(c))
        assert c.keyword == ''
        assert c.value == '       / EXPOSURE INFORMATION'

    def test_specify_undefined_value(self):
        # this is how to specify an undefined value
        c = fits.Card("undef", fits.card.UNDEFINED)
        assert str(c) == _pad("UNDEF   =")

    def test_complex_number_using_string_input(self):
        # complex number using string input
        c = fits.Card.fromstring('abc     = (8, 9)')
        assert str(c) == _pad("ABC     =               (8, 9)")

    def test_fixable_non_standard_fits_card(self):
        # fixable non-standard FITS card will keep the original format
        c = fits.Card.fromstring('abc     = +  2.1   e + 12')
        assert c.value == 2100000000000.0
        assert str(c) == _pad("ABC     =             +2.1E+12")

    def test_fixable_non_fsc(self):
        # fixable non-FSC: if the card is not parsable, it's value will be
        # assumed
        # to be a string and everything after the first slash will be comment
        c = fits.Card.fromstring(
            "no_quote=  this card's value has no quotes "
            "/ let's also try the comment")
        assert (str(c) == "NO_QUOTE= 'this card''s value has no quotes' "
                          "/ let's also try the comment       ")

    def test_undefined_value_using_string_input(self):
        # undefined value using string input
        c = fits.Card.fromstring('abc     =    ')
        assert str(c) == _pad("ABC     =")

    def test_misalocated_equal_sign(self):
        # test mislocated "=" sign
        c = fits.Card.fromstring('xyz= 100')
        assert c.keyword == 'XYZ'
        assert c.value == 100
        assert str(c) == _pad("XYZ     =                  100")

    def test_equal_only_up_to_column_10(self):
        # the test of "=" location is only up to column 10
        c = fits.Card.fromstring("histo       =   (1, 2)")
        assert str(c) == _pad("HISTO   = '=   (1, 2)'")
        c = fits.Card.fromstring("   history          (1, 2)")
        assert str(c) == _pad("HISTO   = 'ry          (1, 2)'")

    def test_verify_invalid_equal_sign(self):
        # verification
        c = fits.Card.fromstring('abc= a6')
        with warnings.catch_warnings(record=True) as w:
            c.verify()
            assert len(w) == 2
            assert ('Card image is not FITS standard (equal sign not at '
                    'column 8)' in str(w[0].message))

    def test_fix_invalid_equal_sign(self):
        c = fits.Card.fromstring('abc= a6')
        with warnings.catch_warnings(record=True) as w:
            c.verify('fix')
            fix_text = 'Fixed card to meet the FITS standard: ABC'
            assert len(w) == 2
            assert fix_text in str(w[0].message)
        assert str(c) == _pad("ABC     = 'a6      '")

    def test_long_string_value(self):
        # test long string value
        c = fits.Card('abc', 'long string value ' * 10, 'long comment ' * 10)
        assert (str(c) ==
            "ABC     = 'long string value long string value long string value long string &' "
            "CONTINUE  'value long string value long string value long string value long &'  "
            "CONTINUE  'string value long string value long string value &'                  "
            "CONTINUE  '&' / long comment long comment long comment long comment long        "
            "CONTINUE  '&' / comment long comment long comment long comment long comment     "
            "CONTINUE  '&' / long comment                                                    ")

    def test_long_string_from_file(self):
        c = fits.Card('abc', 'long string value ' * 10, 'long comment ' * 10)
        hdu = fits.PrimaryHDU()
        hdu.header.append(c)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        c = hdul[0].header.cards['abc']
        hdul.close()
        assert (str(c) ==
            "ABC     = 'long string value long string value long string value long string &' "
            "CONTINUE  'value long string value long string value long string value long &'  "
            "CONTINUE  'string value long string value long string value &'                  "
            "CONTINUE  '&' / long comment long comment long comment long comment long        "
            "CONTINUE  '&' / comment long comment long comment long comment long comment     "
            "CONTINUE  '&' / long comment                                                    ")

    def test_word_in_long_string_too_long(self):
        # if a word in a long string is too long, it will be cut in the middle
        c = fits.Card('abc', 'longstringvalue' * 10, 'longcomment' * 10)
        assert (str(c) ==
            "ABC     = 'longstringvaluelongstringvaluelongstringvaluelongstringvaluelongstr&'"
            "CONTINUE  'ingvaluelongstringvaluelongstringvaluelongstringvaluelongstringvalu&'"
            "CONTINUE  'elongstringvalue&'                                                   "
            "CONTINUE  '&' / longcommentlongcommentlongcommentlongcommentlongcommentlongcomme"
            "CONTINUE  '&' / ntlongcommentlongcommentlongcommentlongcomment                  ")

    def test_long_string_value_via_fromstring(self):
        # long string value via fromstring() method
        c = fits.Card.fromstring(
            _pad("abc     = 'longstring''s testing  &  ' "
                 "/ comments in line 1") +
            _pad("continue  'continue with long string but without the "
                 "ampersand at the end' /") +
            _pad("continue  'continue must have string value (with quotes)' "
                 "/ comments with ''. "))
        assert (str(c) ==
            "ABC     = 'longstring''s testing  continue with long string but without the &'  "
            "CONTINUE  'ampersand at the endcontinue must have string value (with quotes)&'  "
            "CONTINUE  '&' / comments in line 1 comments with ''.                            ")

    def test_continue_card_with_equals_in_value(self):
        """
        Regression test for #117.
        """

        c = fits.Card.fromstring(
            fits.card._pad("EXPR    = '/grp/hst/cdbs//grid/pickles/dat_uvk/pickles_uk_10.fits * &'") +
            fits.card._pad("CONTINUE  '5.87359e-12 * MWAvg(Av=0.12)&'") +
            fits.card._pad("CONTINUE  '&' / pysyn expression"))

        assert c.keyword == 'EXPR'
        assert (c.value ==
                '/grp/hst/cdbs//grid/pickles/dat_uvk/pickles_uk_10.fits '
                '* 5.87359e-12 * MWAvg(Av=0.12)')
        assert c.comment == 'pysyn expression'

    def test_hierarch_card_creation(self):
        # Test automatic upgrade to hierarch card
        with warnings.catch_warnings(record=True) as w:
            c = fits.Card('ESO INS SLIT2 Y1FRML',
                          'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)')
            assert len(w) == 1
            assert 'HIERARCH card will be created' in str(w[0].message)
            assert (str(c) ==
                    "HIERARCH ESO INS SLIT2 Y1FRML= "
                    "'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)'")

        # Test manual creation of hierarch card
        c = fits.Card('hierarch abcdefghi', 10)
        assert (str(c) ==
            "HIERARCH abcdefghi = "
            "10                                                         ")
        c = fits.Card('HIERARCH ESO INS SLIT2 Y1FRML',
                      'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)')
        assert (str(c) ==
            "HIERARCH ESO INS SLIT2 Y1FRML= "
            "'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)'")

    def test_hierarch_card_lookup(self):
        header = fits.Header()
        header['hierarch abcdefghi'] = 10
        assert 'abcdefghi' in header
        assert header['abcdefghi'] == 10
        assert 'ABCDEFGHI' not in header

    def test_header_setitem_invalid(self):
        header = fits.Header()

        def test():
            header['FOO'] = ('bar', 'baz', 'qux')

        pytest.raises(ValueError, test)

    def test_header_setitem_1tuple(self):
        header = fits.Header()
        header['FOO'] = ('BAR',)
        assert header['FOO'] == 'BAR'
        assert header[0] == 'BAR'
        assert header.comments[0] == ''
        assert header.comments['FOO'] == ''

    def test_header_setitem_2tuple(self):
        header = fits.Header()
        header['FOO'] = ('BAR', 'BAZ')
        assert header['FOO'] == 'BAR'
        assert header[0] == 'BAR'
        assert header.comments[0] == 'BAZ'
        assert header.comments['FOO'] == 'BAZ'

    def test_header_set_value_to_none(self):
        """
        Setting the value of a card to None should simply give that card a
        blank value.
        """

        header = fits.Header()
        header['FOO'] = 'BAR'
        assert header['FOO'] == 'BAR'
        header['FOO'] = None
        assert header['FOO'] == ''

    def test_set_comment_only(self):
        header = fits.Header([('A', 'B', 'C')])
        header.set('A', comment='D')
        assert header['A'] == 'B'
        assert header.comments['A'] == 'D'

    def test_header_iter(self):
        header = fits.Header([('A', 'B'), ('C', 'D')])
        assert list(header) == ['A', 'C']

    def test_header_slice(self):
        header = fits.Header([('A', 'B'), ('C', 'D'), ('E', 'F')])
        newheader = header[1:]
        assert len(newheader) == 2
        assert 'A' not in newheader
        assert 'C' in newheader
        assert 'E' in newheader

        newheader = header[::-1]
        assert len(newheader) == 3
        assert newheader[0] == 'F'
        assert newheader[1] == 'D'
        assert newheader[2] == 'B'

        newheader = header[::2]
        assert len(newheader) == 2
        assert 'A' in newheader
        assert 'C' not in newheader
        assert 'E' in newheader

    def test_header_slice_assignment(self):
        """
        Assigning to a slice should just assign new values to the cards
        included in the slice.
        """

        header = fits.Header([('A', 'B'), ('C', 'D'), ('E', 'F')])

        # Test assigning slice to the same value; this works similarly to numpy
        # arrays
        header[1:] = 1
        assert header[1] == 1
        assert header[2] == 1

        # Though strings are iterable they should be treated as a scalar value
        header[1:] = 'GH'
        assert header[1] == 'GH'
        assert header[2] == 'GH'

        # Now assign via an iterable
        header[1:] = ['H', 'I']
        assert header[1] == 'H'
        assert header[2] == 'I'

    def test_header_slice_delete(self):
        """Test deleting a slice of cards from the header."""

        header = fits.Header([('A', 'B'), ('C', 'D'), ('E', 'F')])
        del header[1:]
        assert len(header) == 1
        assert header[0] == 'B'
        del header[:]
        assert len(header) == 0

    def test_wildcard_slice(self):
        """Test selecting a subsection of a header via wildcard matching."""

        header = fits.Header([('ABC', 0), ('DEF', 1), ('ABD', 2)])
        newheader = header['AB*']
        assert len(newheader) == 2
        assert newheader[0] == 0
        assert newheader[1] == 2

    def test_wildcard_slice_assignment(self):
        """Test assigning to a header slice selected via wildcard matching."""

        header = fits.Header([('ABC', 0), ('DEF', 1), ('ABD', 2)])

        # Test assigning slice to the same value; this works similarly to numpy
        # arrays
        header['AB*'] = 1
        assert header[0] == 1
        assert header[2] == 1

        # Though strings are iterable they should be treated as a scalar value
        header['AB*'] = 'GH'
        assert header[0] == 'GH'
        assert header[2] == 'GH'

        # Now assign via an iterable
        header['AB*'] = ['H', 'I']
        assert header[0] == 'H'
        assert header[2] == 'I'

    def test_wildcard_slice_deletion(self):
        """Test deleting cards from a header that match a wildcard pattern."""

        header = fits.Header([('ABC', 0), ('DEF', 1), ('ABD', 2)])
        del header['AB*']
        assert len(header) == 1
        assert header[0] == 1

    def test_header_history(self):
        header = fits.Header([('ABC', 0), ('HISTORY', 1), ('HISTORY', 2),
                              ('DEF', 3), ('HISTORY', 4), ('HISTORY', 5)])
        assert header['HISTORY'] == [1, 2, 4, 5]

    def test_header_clear(self):
        header = fits.Header([('A', 'B'), ('C', 'D')])
        header.clear()
        assert 'A' not in header
        assert 'C' not in header
        assert len(header) == 0

    def test_header_fromkeys(self):
        header = fits.Header.fromkeys(['A', 'B'])
        assert 'A' in header
        assert header['A'] == ''
        assert header.comments['A'] == ''
        assert 'B' in header
        assert header['B'] == ''
        assert header.comments['B'] == ''

    def test_header_fromkeys_with_value(self):
        header = fits.Header.fromkeys(['A', 'B'], 'C')
        assert 'A' in header
        assert header['A'] == 'C'
        assert header.comments['A'] == ''
        assert 'B' in header
        assert header['B'] == 'C'
        assert header.comments['B'] == ''

    def test_header_fromkeys_with_value_and_comment(self):
        header = fits.Header.fromkeys(['A'], ('B', 'C'))
        assert 'A' in header
        assert header['A'] == 'B'
        assert header.comments['A'] == 'C'

    def test_header_fromkeys_with_duplicates(self):
        header = fits.Header.fromkeys(['A', 'B', 'A'], 'C')
        assert 'A' in header
        assert ('A', 0) in header
        assert ('A', 1) in header
        assert ('A', 2) not in header
        assert header[0] == 'C'
        assert header['A'] == 'C'
        assert header[('A', 0)] == 'C'
        assert header[2] == 'C'
        assert header[('A', 1)] == 'C'

    def test_header_items(self):
        header = fits.Header([('A', 'B'), ('C', 'D')])
        assert header.items() == list(header.iteritems())

    def test_header_iterkeys(self):
        header = fits.Header([('A', 'B'), ('C', 'D')])
        for a, b in itertools.izip(header.iterkeys(), header):
            assert a == b

    def test_header_itervalues(self):
        header = fits.Header([('A', 'B'), ('C', 'D')])
        for a, b in itertools.izip(header.itervalues(), ['B', 'D']):
            assert a == b

    def test_header_keys(self):
        hdul = fits.open(self.data('arange.fits'))
        assert (hdul[0].header.keys() ==
                ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NAXIS3',
                 'EXTEND'])

    def test_header_list_like_pop(self):
        header = fits.Header([('A', 'B'), ('C', 'D'), ('E', 'F'),
                              ('G', 'H')])

        last = header.pop()
        assert last == 'H'
        assert len(header) == 3
        assert header.keys() == ['A', 'C', 'E']

        mid = header.pop(1)
        assert mid == 'D'
        assert len(header) == 2
        assert header.keys() == ['A', 'E']

        first = header.pop(0)
        assert first == 'B'
        assert len(header) == 1
        assert header.keys() == ['E']

        pytest.raises(IndexError, header.pop, 42)

    def test_header_dict_like_pop(self):
        header = fits.Header([('A', 'B'), ('C', 'D'), ('E', 'F'),
                              ('G', 'H')])
        pytest.raises(TypeError, header.pop, 'A', 'B', 'C')

        last = header.pop('G')
        assert last == 'H'
        assert len(header) == 3
        assert header.keys() == ['A', 'C', 'E']

        mid = header.pop('C')
        assert mid == 'D'
        assert len(header) == 2
        assert header.keys() == ['A', 'E']

        first = header.pop('A')
        assert first == 'B'
        assert len(header) == 1
        assert header.keys() == ['E']

        default = header.pop('X', 'Y')
        assert default == 'Y'
        assert len(header) == 1

        pytest.raises(KeyError, header.pop, 'X')

    def test_popitem(self):
        header = fits.Header([('A', 'B'), ('C', 'D'), ('E', 'F')])
        keyword, value = header.popitem()
        assert keyword not in header
        assert len(header) == 2
        keyword, value = header.popitem()
        assert keyword not in header
        assert len(header) == 1
        keyword, value = header.popitem()
        assert keyword not in header
        assert len(header) == 0
        pytest.raises(KeyError, header.popitem)

    def test_setdefault(self):
        header = fits.Header([('A', 'B'), ('C', 'D'), ('E', 'F')])
        assert header.setdefault('A') == 'B'
        assert header.setdefault('C') == 'D'
        assert header.setdefault('E') == 'F'
        assert len(header) == 3
        assert header.setdefault('G', 'H') == 'H'
        assert len(header) == 4
        assert 'G' in header
        assert header.setdefault('G', 'H') == 'H'
        assert len(header) == 4

    def test_update_from_dict(self):
        """
        Test adding new cards and updating existing cards from a dict using
        Header.update()
        """

        header = fits.Header([('A', 'B'), ('C', 'D')])
        header.update({'A': 'E', 'F': 'G'})
        assert header['A'] == 'E'
        assert header[0] == 'E'
        assert 'F' in header
        assert header['F'] == 'G'
        assert header[-1] == 'G'

        # Same as above but this time pass the update dict as keyword arguments
        header = fits.Header([('A', 'B'), ('C', 'D')])
        header.update(A='E', F='G')
        assert header['A'] == 'E'
        assert header[0] == 'E'
        assert 'F' in header
        assert header['F'] == 'G'
        assert header[-1] == 'G'

    def test_update_from_iterable(self):
        """
        Test adding new cards and updating existing cards from an iterable of
        cards and card tuples.
        """

        header = fits.Header([('A', 'B'), ('C', 'D')])
        header.update([('A', 'E'), fits.Card('F', 'G')])
        assert header['A'] == 'E'
        assert header[0] == 'E'
        assert 'F' in header
        assert header['F'] == 'G'
        assert header[-1] == 'G'

    def test_header_extend(self):
        """
        Test extending a header both with and without stripping cards from the
        extension header.
        """

        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu2.header['MYKEY'] = ('some val', 'some comment')
        hdu.header += hdu2.header
        assert len(hdu.header) == 5
        assert hdu.header[-1] == 'some val'

        # Same thing, but using + instead of +=
        hdu = fits.PrimaryHDU()
        hdu.header = hdu.header + hdu2.header
        assert len(hdu.header) == 5
        assert hdu.header[-1] == 'some val'

        # Directly append the other header in full--not usually a desirable
        # operation when the header is coming from another HDU
        hdu.header.extend(hdu2.header, strip=False)
        assert len(hdu.header) == 11
        assert hdu.header.keys()[5] == 'XTENSION'
        assert hdu.header[-1] == 'some val'
        assert ('MYKEY', 1) in hdu.header

    def test_header_extend_unique(self):
        """
        Test extending the header with and without unique=True.
        """
        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu.header['MYKEY'] = ('some val', 'some comment')
        hdu2.header['MYKEY'] = ('some other val', 'some other comment')
        hdu.header.extend(hdu2.header)
        assert len(hdu.header) == 6
        assert hdu.header[-2] == 'some val'
        assert hdu.header[-1] == 'some other val'

        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu.header['MYKEY'] = ('some val', 'some comment')
        hdu.header.extend(hdu2.header, unique=True)
        assert len(hdu.header) == 5
        assert hdu.header[-1] == 'some val'

    def test_header_extend_unique(self):
        """
        Test extending the header with and without unique=True.
        """
        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu.header['MYKEY'] = ('some val', 'some comment')
        hdu2.header['MYKEY'] = ('some other val', 'some other comment')
        hdu.header.extend(hdu2.header)
        assert len(hdu.header) == 6
        assert hdu.header[-2] == 'some val'
        assert hdu.header[-1] == 'some other val'

        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu.header['MYKEY'] = ('some val', 'some comment')
        hdu.header.extend(hdu2.header, unique=True)
        assert len(hdu.header) == 5
        assert hdu.header[-1] == 'some val'

    def test_header_extend_update(self):
        """
        Test extending the header with and without update=True.
        """

        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu.header['MYKEY'] = ('some val', 'some comment')
        hdu.header['HISTORY'] = 'history 1'
        hdu2.header['MYKEY'] = ('some other val', 'some other comment')
        hdu2.header['HISTORY'] = 'history 1'
        hdu2.header['HISTORY'] = 'history 2'
        hdu.header.extend(hdu2.header)
        assert len(hdu.header) == 9
        assert ('MYKEY', 0) in hdu.header
        assert ('MYKEY', 1) in hdu.header
        assert hdu.header[('MYKEY', 1)] == 'some other val'
        assert len(hdu.header['HISTORY']) == 3
        assert hdu.header[-1] == 'history 2'

        hdu = fits.PrimaryHDU()
        hdu.header['MYKEY'] = ('some val', 'some comment')
        hdu.header['HISTORY'] = 'history 1'
        hdu.header.extend(hdu2.header, update=True)
        assert len(hdu.header) == 7
        assert ('MYKEY', 0) in hdu.header
        assert ('MYKEY', 1) not in hdu.header
        assert hdu.header['MYKEY'] == 'some other val'
        assert len(hdu.header['HISTORY']) == 2
        assert hdu.header[-1] == 'history 2'

    def test_header_count(self):
        header = fits.Header([('A', 'B'), ('C', 'D'), ('E', 'F')])
        assert header.count('A') == 1
        assert header.count('C') == 1
        assert header.count('E') == 1
        header['HISTORY'] = 'a'
        header['HISTORY'] = 'b'
        assert header.count('HISTORY') == 2
        pytest.raises(KeyError, header.count, 'G')

    def test_header_append_use_blanks(self):
        """
        Tests that blank cards can be appended, and that future appends will
        use blank cards when available (unless useblanks=False)
        """

        header = fits.Header([('A', 'B'), ('C', 'D')])

        # Append a couple blanks
        header.append()
        header.append()
        assert len(header) == 4
        assert header[-1] == ''
        assert header[-2] == ''

        # New card should fill the first blank by default
        header.append(('E', 'F'))
        assert len(header) == 4
        assert header[-2] == 'F'
        assert header[-1] == ''

        # This card should not use up a blank spot
        header.append(('G', 'H'), useblanks=False)
        assert len(header) == 5
        assert header[-1] == ''
        assert header[-2] == 'H'

    def test_header_append_keyword_only(self):
        """
        Test appending a new card with just the keyword, and no value or
        comment given.
        """

        header = fits.Header([('A', 'B'), ('C', 'D')])
        header.append('E')
        assert len(header) == 3
        assert header.keys()[-1] == 'E'
        assert header[-1] == ''
        assert header.comments['E'] == ''

        # Try appending a blank--normally this can be accomplished with just
        # header.append(), but header.append('') should also work (and is maybe
        # a little more clear)
        header.append('')
        assert len(header) == 4
        # Blank keywords are ignored in the keys list
        assert header.keys()[-1] == 'E'
        assert header[''] == ''

    def test_header_insert_use_blanks(self):
        header = fits.Header([('A', 'B'), ('C', 'D')])

        # Append a couple blanks
        header.append()
        header.append()

        # Insert a new card; should use up one of the blanks
        header.insert(1, ('E', 'F'))
        assert len(header) == 4
        assert header[1] == 'F'
        assert header[-1] == ''
        assert header[-2] == 'D'

        # Insert a new card without using blanks
        header.insert(1, ('G', 'H'), useblanks=False)
        assert len(header) == 5
        assert header[1] == 'H'
        assert header[-1] == ''

    def test_header_comments(self):
        header = fits.Header([('A', 'B', 'C'), ('DEF', 'G', 'H')])
        assert (repr(header.comments) ==
                '       A  C\n'
                '     DEF  H')

    def test_comment_slices_and_filters(self):
        header = fits.Header([('AB', 'C', 'D'), ('EF', 'G', 'H'),
                              ('AI', 'J', 'K')])
        s = header.comments[1:]
        assert list(s) == ['H', 'K']
        s = header.comments[::-1]
        assert list(s) == ['K', 'H', 'D']
        s = header.comments['A*']
        assert list(s) == ['D', 'K']

    def test_comment_slice_filter_assign(self):
        header = fits.Header([('AB', 'C', 'D'), ('EF', 'G', 'H'),
                              ('AI', 'J', 'K')])
        header.comments[1:] = 'L'
        assert list(header.comments) == ['D', 'L', 'L']
        assert header.cards[header.index('AB')].comment == 'D'
        assert header.cards[header.index('EF')].comment == 'L'
        assert header.cards[header.index('AI')].comment == 'L'

        header.comments[::-1] = header.comments[:]
        assert list(header.comments) == ['L', 'L', 'D']

        header.comments['A*'] = ['M', 'N']
        assert list(header.comments) == ['M', 'L', 'N']

    def test_update_comment(self):
        hdul = fits.open(self.data('arange.fits'))
        hdul[0].header['FOO'] = ('BAR', 'BAZ')
        hdul.writeto(self.temp('test.fits'))

        hdul = fits.open(self.temp('test.fits'), mode='update')
        hdul[0].header.comments['FOO'] = 'QUX'
        hdul.close()

        hdul = fits.open(self.temp('test.fits'))
        assert hdul[0].header.comments['FOO'] == 'QUX'

    def test_update_commentary(self):
        header = fits.Header()
        header['FOO'] = 'BAR'
        header['HISTORY'] = 'ABC'
        header['FRED'] = 'BARNEY'
        header['HISTORY'] = 'DEF'
        header['HISTORY'] = 'GHI'

        assert header['HISTORY'], ['ABC', 'DEF' == 'GHI']

        # Single value update
        header['HISTORY'][0] = 'FOO'
        assert header['HISTORY'], ['FOO', 'DEF' == 'GHI']

        # Single value partial slice update
        header['HISTORY'][1:] = 'BAR'
        assert header['HISTORY'], ['FOO', 'BAR' == 'BAR']

        # Multi-value update
        header['HISTORY'][:] = ['BAZ', 'QUX']
        assert header['HISTORY'], ['BAZ', 'QUX' == 'BAR']

    def test_long_commentary_card(self):
        header = fits.Header()
        header['FOO'] = 'BAR'
        header['BAZ'] = 'QUX'
        longval = 'ABC' * 30
        header['HISTORY'] = longval
        header['FRED'] = 'BARNEY'
        header['HISTORY'] = longval

        assert len(header) == 7
        assert header.keys()[2] == 'FRED'
        assert str(header.cards[3]) == 'HISTORY ' + longval[:72]
        assert str(header.cards[4]).rstrip() == 'HISTORY ' + longval[72:]

        header.set('HISTORY', longval, after='FOO')
        assert len(header) == 9
        assert str(header.cards[1]) == 'HISTORY ' + longval[:72]
        assert str(header.cards[2]).rstrip() == 'HISTORY ' + longval[72:]

    def test_header_fromtextfile(self):
        """Regression test for #122.

        Manually write a text file containing some header cards ending with
        newlines and ensure that fromtextfile can read them back in.
        """

        header = fits.Header()
        header['A'] = ('B', 'C')
        header['B'] = ('C', 'D')
        header['C'] = ('D', 'E')

        with open(self.temp('test.hdr'), 'w') as f:
            f.write('\n'.join(str(c).strip() for c in header.cards))

        header2 = fits.Header.fromtextfile(self.temp('test.hdr'))
        assert header == header2


class TestRecordValuedKeywordCards(FitsTestCase):
    """
    Tests for handling of record-valued keyword cards as used by the FITS WCS
    Paper IV proposal.

    These tests are derived primarily from the release notes for PyFITS 1.4 (in
    which this feature was first introduced.
    """

    def setup(self):
        super(TestRecordValuedKeywordCards, self).setup()
        self._test_header = fits.Header()
        self._test_header.set('DP1', 'NAXIS: 2')
        self._test_header.set('DP1', 'AXIS.1: 1')
        self._test_header.set('DP1', 'AXIS.2: 2')
        self._test_header.set('DP1', 'NAUX: 2')
        self._test_header.set('DP1', 'AUX.1.COEFF.0: 0')
        self._test_header.set('DP1', 'AUX.1.POWER.0: 1')
        self._test_header.set('DP1', 'AUX.1.COEFF.1: 0.00048828125')
        self._test_header.set('DP1', 'AUX.1.POWER.1: 1')

    def test_initialize_rvkc(self):
        """
        Test different methods for initializing a card that should be
        recognized as a RVKC
        """

        c = fits.Card.fromstring("DP1     = 'NAXIS: 2' / A comment")
        assert c.keyword == 'DP1.NAXIS'
        assert c.value == 2.0
        assert c.field_specifier == 'NAXIS'
        assert c.comment == 'A comment'

        c = fits.Card.fromstring("DP1     = 'NAXIS: 2.1'")
        assert c.keyword == 'DP1.NAXIS'
        assert c.value == 2.1
        assert c.field_specifier == 'NAXIS'

        c = fits.Card.fromstring("DP1     = 'NAXIS: a'")
        assert c.keyword == 'DP1'
        assert c.value == 'NAXIS: a'
        assert c.field_specifier == None

        c = fits.Card('DP1', 'NAXIS: 2')
        assert c.keyword == 'DP1.NAXIS'
        assert c.value == 2.0
        assert c.field_specifier == 'NAXIS'

        c = fits.Card('DP1', 'NAXIS: 2.0')
        assert c.keyword == 'DP1.NAXIS'
        assert c.value == 2.0
        assert c.field_specifier == 'NAXIS'

        c = fits.Card('DP1', 'NAXIS: a')
        assert c.keyword == 'DP1'
        assert c.value == 'NAXIS: a'
        assert c.field_specifier == None

        c = fits.Card('DP1.NAXIS', 2)
        assert c.keyword == 'DP1.NAXIS'
        assert c.value == 2.0
        assert c.field_specifier == 'NAXIS'

        c = fits.Card('DP1.NAXIS', 2.0)
        assert c.keyword == 'DP1.NAXIS'
        assert c.value == 2.0
        assert c.field_specifier == 'NAXIS'

        with ignore_warnings():
            c = fits.Card('DP1.NAXIS', 'a')
        assert c.keyword == 'DP1.NAXIS'
        assert c.value == 'a'
        assert c.field_specifier == None

    def test_parse_field_specifier(self):
        """
        Tests that the field_specifier can accessed from a card read from a
        string before any other attributes are accessed.
        """

        c = fits.Card.fromstring("DP1     = 'NAXIS: 2' / A comment")
        assert c.field_specifier == 'NAXIS'
        assert c.keyword == 'DP1.NAXIS'
        assert c.value == 2.0
        assert c.comment == 'A comment'

    def test_update_field_specifier(self):
        """
        Test setting the field_specifier attribute and updating the card image
        to reflect the new value.
        """

        c = fits.Card.fromstring("DP1     = 'NAXIS: 2' / A comment")
        assert c.field_specifier == 'NAXIS'
        c.field_specifier = 'NAXIS1'
        assert c.field_specifier == 'NAXIS1'
        assert c.keyword == 'DP1.NAXIS1'
        assert c.value == 2.0
        assert c.comment == 'A comment'
        assert str(c).rstrip() == "DP1     = 'NAXIS1: 2' / A comment"

    def test_field_specifier_case_senstivity(self):
        """
        The keyword portion of an RVKC should still be case-insensitive, but
        the field-specifier portion should be case-sensitive.
        """

        header = fits.Header()
        header.set('abc.def', 1)
        header.set('abc.DEF', 2)
        assert header['abc.def'] == 1
        assert header['ABC.def'] == 1
        assert header['aBc.def'] == 1
        assert header['ABC.DEF'] == 2
        assert 'ABC.dEf' not in header

    def test_get_rvkc_by_index(self):
        """
        Returning a RVKC from a header via index lookup should return the
        float value of the card.
        """

        assert self._test_header[0] == 2.0
        assert isinstance(self._test_header[0], float)
        assert self._test_header[1] == 1.0
        assert isinstance(self._test_header[1], float)

    def test_get_rvkc_by_keyword(self):
        """
        Returning a RVKC just via the keyword name should return the floating
        point value of the first card with that keyword.
        """

        assert self._test_header['DP1'] == 2.0

    def test_get_rvkc_by_keyword_and_field_specifier(self):
        """
        Returning a RVKC via the full keyword/field-specifier combination
        should return the floating point value associated with the RVKC.
        """

        assert self._test_header['DP1.NAXIS'] == 2.0
        assert isinstance(self._test_header['DP1.NAXIS'], float)
        assert self._test_header['DP1.AUX.1.COEFF.1'] == 0.00048828125

    def test_access_nonexistent_rvkc(self):
        """
        Accessing a nonexistent RVKC should raise an IndexError for
        index-based lookup, or a KeyError for keyword lookup (like a normal
        card).
        """

        pytest.raises(IndexError, lambda x: self._test_header[x], 8)
        pytest.raises(KeyError, lambda k: self._test_header[k], 'DP1.AXIS.3')

    def test_update_rvkc(self):
        """A RVKC can be updated either via index or keyword access."""

        self._test_header[0] = 3
        assert self._test_header['DP1.NAXIS'] == 3.0
        assert isinstance(self._test_header['DP1.NAXIS'], float)

        self._test_header['DP1.AXIS.1'] = 1.1
        assert self._test_header['DP1.AXIS.1'] == 1.1

    def test_rvkc_insert_after(self):
        """
        It should be possible to insert a new RVKC after an existing one
        specified by the full keyword/field-specifier combination."""

        self._test_header.set('DP1', 'AXIS.3: 1', 'a comment',
                              after='DP1.AXIS.2')
        assert self._test_header[3] == 1
        assert self._test_header['DP1.AXIS.3'] == 1

    def test_rvkc_delete(self):
        """
        Deleting a RVKC should work as with a normal card by using the full
        keyword/field-spcifier combination.
        """

        del self._test_header['DP1.AXIS.1']
        assert len(self._test_header) == 7
        assert self._test_header.keys()[0] == 'DP1.NAXIS'
        assert self._test_header[0] == 2
        assert self._test_header.keys()[1] == 'DP1.AXIS.2'
        assert self._test_header[1] == 2

    def test_pattern_matching_keys(self):
        """Test the keyword filter strings with RVKCs."""

        cl = self._test_header['DP1.AXIS.*']
        assert isinstance(cl, fits.Header)
        assert ([str(c).strip() for c in cl.cards] ==
                ["DP1     = 'AXIS.1: 1'",
                 "DP1     = 'AXIS.2: 2'"])

        cl = self._test_header['DP1.N*']
        assert ([str(c).strip() for c in cl.cards] ==
                ["DP1     = 'NAXIS: 2'",
                 "DP1     = 'NAUX: 2'"])

        cl = self._test_header['DP1.AUX...']
        assert ([str(c).strip() for c in cl.cards] ==
                ["DP1     = 'AUX.1.COEFF.0: 0'",
                 "DP1     = 'AUX.1.POWER.0: 1'",
                 "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
                 "DP1     = 'AUX.1.POWER.1: 1'"])

        cl = self._test_header['DP?.NAXIS']
        assert ([str(c).strip() for c in cl.cards] ==
                ["DP1     = 'NAXIS: 2'"])

        cl = self._test_header['DP1.A*S.*']
        assert ([str(c).strip() for c in cl.cards] ==
                ["DP1     = 'AXIS.1: 1'",
                 "DP1     = 'AXIS.2: 2'"])

    def test_pattern_matching_key_deletion(self):
        """Deletion by filter strings should work."""

        del self._test_header['DP1.A*...']
        assert len(self._test_header) == 2
        assert self._test_header.keys()[0] == 'DP1.NAXIS'
        assert self._test_header[0] == 2
        assert self._test_header.keys()[1] == 'DP1.NAUX'
        assert self._test_header[1] == 2

    def test_successive_pattern_matching(self):
        """
        A card list returned via a filter string should be further filterable.
        """

        cl = self._test_header['DP1.A*...']
        assert ([str(c).strip() for c in cl.cards] ==
                ["DP1     = 'AXIS.1: 1'",
                 "DP1     = 'AXIS.2: 2'",
                 "DP1     = 'AUX.1.COEFF.0: 0'",
                 "DP1     = 'AUX.1.POWER.0: 1'",
                 "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
                 "DP1     = 'AUX.1.POWER.1: 1'"])

        cl2 = cl['*.*AUX...']
        assert ([str(c).strip() for c in cl2.cards] ==
                ["DP1     = 'AUX.1.COEFF.0: 0'",
                 "DP1     = 'AUX.1.POWER.0: 1'",
                 "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
                 "DP1     = 'AUX.1.POWER.1: 1'"])

    def test_rvkc_in_cardlist_keys(self):
        """
        The CardList.keys() method should return full keyword/field-spec values
        for RVKCs.
        """

        cl = self._test_header['DP1.AXIS.*']
        assert cl.keys() == ['DP1.AXIS.1', 'DP1.AXIS.2']

    def test_rvkc_in_cardlist_values(self):
        """
        The CardList.values() method should return the values of all RVKCs as
        floating point values.
        """

        cl = self._test_header['DP1.AXIS.*']
        assert cl.values() == [1.0, 2.0]

    def test_rvkc_value_attribute(self):
        """
        Individual card values should be accessible by the .value attribute
        (which should return a float).
        """

        cl = self._test_header['DP1.AXIS.*']
        assert cl.cards[0].value == 1.0
        assert isinstance(cl.cards[0].value, float)
