# Licensed under a 3-clause BSD style license - see PYFITS.rst

import collections
import copy
import warnings
from io import BytesIO, StringIO

import numpy as np
import pytest

from astropy.io import fits
from astropy.io.fits.card import _pad
from astropy.io.fits.header import _pad_length
from astropy.io.fits.util import encode_ascii
from astropy.io.fits.verify import VerifyError, VerifyWarning
from astropy.utils.exceptions import AstropyUserWarning
from astropy.utils.misc import _NOT_OVERWRITING_MSG_MATCH

from .conftest import FitsTestCase


def test_shallow_copy():
    """Make sure that operations on a shallow copy do not alter the original.
    #4990."""
    original_header = fits.Header([("a", 1), ("b", 1)])
    copied_header = copy.copy(original_header)

    # Modifying the original dict should not alter the copy
    original_header["c"] = 100
    assert "c" not in copied_header

    # and changing the copy should not change the original.
    copied_header["a"] = 0
    assert original_header["a"] == 1


def test_init_with_header():
    """Make sure that creating a Header from another Header makes a copy if
    copy is True."""

    original_header = fits.Header([("a", 10)])
    new_header = fits.Header(original_header, copy=True)
    original_header["a"] = 20
    assert new_header["a"] == 10

    new_header["a"] = 0
    assert original_header["a"] == 20


def test_init_with_dict():
    dict1 = {"a": 11, "b": 12, "c": 13, "d": 14, "e": 15}
    h1 = fits.Header(dict1)
    for i in dict1:
        assert dict1[i] == h1[i]


def test_init_with_ordereddict():
    # Create a list of tuples. Each tuple consisting of a letter and the number
    list1 = [(i, j) for j, i in enumerate("abcdefghijklmnopqrstuvwxyz")]
    # Create an ordered dictionary and a header from this dictionary
    dict1 = collections.OrderedDict(list1)
    h1 = fits.Header(dict1)
    # Check that the order is preserved of the initial list
    assert all(h1[val] == list1[i][1] for i, val in enumerate(h1))


class TestHeaderFunctions(FitsTestCase):
    """Test Header and Card objects."""

    def test_rename_keyword(self):
        """Test renaming keyword with rename_keyword."""
        header = fits.Header([("A", "B", "C"), ("D", "E", "F")])
        header.rename_keyword("A", "B")
        assert "A" not in header
        assert "B" in header
        assert header[0] == "B"
        assert header["B"] == "B"
        assert header.comments["B"] == "C"

    @pytest.mark.parametrize("key", ["A", "a"])
    def test_indexing_case(self, key):
        """Check that indexing is case insensitive"""
        header = fits.Header([("A", "B", "C"), ("D", "E", "F")])
        assert key in header
        assert header[key] == "B"
        assert header.get(key) == "B"
        assert header.index(key) == 0
        assert header.comments[key] == "C"
        assert header.count(key) == 1
        header.remove(key, ignore_missing=False)

    def test_card_constructor_default_args(self):
        """Test Card constructor with default argument values."""

        c = fits.Card()
        assert c.keyword == ""

    def test_card_from_bytes(self):
        """
        Test loading a Card from a `bytes` object (assuming latin-1 encoding).
        """

        c = fits.Card.fromstring(b"ABC     = 'abc'")
        assert c.keyword == "ABC"
        assert c.value == "abc"

    def test_string_value_card(self):
        """Test Card constructor with string value"""

        c = fits.Card("abc", "<8 ch")
        assert str(c) == _pad("ABC     = '<8 ch   '")
        c = fits.Card("nullstr", "")
        assert str(c) == _pad("NULLSTR = ''")

    def test_boolean_value_card(self):
        """Test Card constructor with boolean value"""

        c = fits.Card("abc", True)
        assert str(c) == _pad("ABC     =                    T")

        c = fits.Card.fromstring("ABC     = F")
        assert c.value is False

    def test_long_integer_value_card(self):
        """Test Card constructor with long integer value"""

        c = fits.Card("long_int", -467374636747637647347374734737437)
        assert str(c) == _pad("LONG_INT= -467374636747637647347374734737437")

    def test_floating_point_value_card(self):
        """Test Card constructor with floating point value"""

        c = fits.Card("floatnum", -467374636747637647347374734737437.0)

        if str(c) != _pad("FLOATNUM= -4.6737463674763E+32") and str(c) != _pad(
            "FLOATNUM= -4.6737463674763E+032"
        ):
            assert str(c) == _pad("FLOATNUM= -4.6737463674763E+32")

    def test_floating_point_string_representation_card(self):
        """
        Ensures Card formats float values with the correct precision, avoiding
        comment truncation

        Regression test for https://github.com/astropy/astropy/issues/14507
        """
        k = "HIERARCH ABC DEF GH IJKLMN"
        com = "[m] abcdef ghijklm nopqrstu vw xyzab"
        c = fits.Card(k, 0.009125, com)
        expected_str = f"{k} = 0.009125 / {com}"
        assert str(c)[: len(expected_str)] == expected_str

        c = fits.Card(k, 8.95, com)
        expected_str = f"{k} = 8.95 / {com}"
        assert str(c)[: len(expected_str)] == expected_str

        c = fits.Card(k, -99.9, com)
        expected_str = f"{k} = -99.9 / {com}"
        assert str(c)[: len(expected_str)] == expected_str

    def test_complex_value_card(self):
        """Test Card constructor with complex value"""

        c = fits.Card("abc", (1.2345377437887837487e88 + 6324767364763746367e-33j))
        f1 = _pad("ABC     = (1.23453774378878E+88, 6.32476736476374E-15)")
        f2 = _pad("ABC     = (1.2345377437887E+088, 6.3247673647637E-015)")
        f3 = _pad("ABC     = (1.23453774378878E+88, 6.32476736476374E-15)")
        if str(c) != f1 and str(c) != f2:
            assert str(c) == f3

    def test_card_image_constructed_too_long(self):
        """Test that over-long cards truncate the comment"""

        # card image constructed from key/value/comment is too long
        # (non-string value)
        c = fits.Card("abc", 9, "abcde" * 20)
        with pytest.warns(fits.verify.VerifyWarning):
            assert (
                str(c) == "ABC     =                    9 "
                "/ abcdeabcdeabcdeabcdeabcdeabcdeabcdeabcdeabcdeab"
            )
        c = fits.Card("abc", "a" * 68, "abcdefg")
        with pytest.warns(fits.verify.VerifyWarning):
            assert str(c) == f"ABC     = '{'a' * 68}'"

    def test_constructor_filter_illegal_data_structures(self):
        """Test that Card constructor raises exceptions on bad arguments"""

        pytest.raises(ValueError, fits.Card, ("abc",), {"value": (2, 3)})
        pytest.raises(ValueError, fits.Card, "key", [], "comment")

    def test_keyword_too_long(self):
        """Test that long Card keywords are allowed, but with a warning"""

        pytest.warns(UserWarning, fits.Card, "abcdefghi", "long")

    def test_illegal_characters_in_key(self):
        """
        Test that Card constructor allows illegal characters in the keyword,
        but creates a HIERARCH card.
        """

        # This test used to check that a ValueError was raised, because a
        # keyword like 'abc+' was simply not allowed.  Now it should create a
        # HIERARCH card.

        with pytest.warns(AstropyUserWarning) as w:
            c = fits.Card("abc+", 9)
        assert len(w) == 1
        assert c.image == _pad("HIERARCH abc+ =                    9")

    def test_add_history(self):
        header = fits.Header(
            [
                ("A", "B", "C"),
                ("HISTORY", 1),
                ("HISTORY", 2),
                ("HISTORY", 3),
                ("", "", ""),
                ("", "", ""),
            ]
        )
        header.add_history(4)
        # One of the blanks should get used, so the length shouldn't change
        assert len(header) == 6
        assert header.cards[4].value == 4
        assert header["HISTORY"] == [1, 2, 3, 4]
        assert repr(header["HISTORY"]) == "1\n2\n3\n4"

        header.add_history(0, after="A")
        assert len(header) == 6
        assert header.cards[1].value == 0
        assert header["HISTORY"] == [0, 1, 2, 3, 4]

    def test_add_blank(self):
        header = fits.Header(
            [("A", "B", "C"), ("", 1), ("", 2), ("", 3), ("", "", ""), ("", "", "")]
        )
        header.add_blank(4)
        # This time a new blank should be added, and the existing blanks don't
        # get used... (though this is really kinda sketchy--there's a
        # distinction between truly blank cards, and cards with blank keywords
        # that isn't currently made int he code)
        assert len(header) == 7
        assert header.cards[6].value == 4
        assert header[""] == [1, 2, 3, "", "", 4]
        assert repr(header[""]) == "1\n2\n3\n\n\n4"

        header.add_blank(0, after="A")
        assert len(header) == 8
        assert header.cards[1].value == 0
        assert header[""] == [0, 1, 2, 3, "", "", 4]

        header[""] = 5
        header[" "] = 6
        assert header[""] == [0, 1, 2, 3, "", "", 4, 5, 6]
        assert header[" "] == [0, 1, 2, 3, "", "", 4, 5, 6]

    def test_update(self):
        class FakeHeader(list):
            def keys(self):
                return [l[0] for l in self]

            def __getitem__(self, key):
                return next(l[1:] for l in self if l[0] == key)

        header = fits.Header()
        header.update({"FOO": ("BAR", "BAZ")})
        header.update(FakeHeader([("A", 1), ("B", 2, "comment")]))

        assert set(header.keys()) == {"FOO", "A", "B"}
        assert header.comments["B"] == "comment"

        # test that comments are preserved
        tmphdr = fits.Header()
        tmphdr["HELLO"] = (1, "this is a comment")
        header.update(tmphdr)
        assert set(header.keys()) == {"FOO", "A", "B", "HELLO"}
        assert header.comments["HELLO"] == "this is a comment"

        header.update(NAXIS1=100, NAXIS2=100)
        assert set(header.keys()) == {"FOO", "A", "B", "HELLO", "NAXIS1", "NAXIS2"}
        assert tuple(header.values()) == ("BAR", 1, 2, 1, 100, 100)

    def test_update_comment(self):
        hdul = fits.open(self.data("arange.fits"))
        hdul[0].header.update({"FOO": ("BAR", "BAZ")})
        assert hdul[0].header["FOO"] == "BAR"
        assert hdul[0].header.comments["FOO"] == "BAZ"

        with pytest.raises(ValueError):
            hdul[0].header.update({"FOO2": ("BAR", "BAZ", "EXTRA")})

        hdul.writeto(self.temp("test.fits"))
        hdul.close()

        hdul = fits.open(self.temp("test.fits"), mode="update")
        hdul[0].header.comments["FOO"] = "QUX"
        hdul.close()

        hdul = fits.open(self.temp("test.fits"))
        assert hdul[0].header.comments["FOO"] == "QUX"

        hdul[0].header.add_comment(0, after="FOO")
        assert str(hdul[0].header.cards[-1]).strip() == "COMMENT 0"
        hdul.close()

    def test_commentary_cards(self):
        # commentary cards
        val = "A commentary card's value has no quotes around it."
        c = fits.Card("HISTORY", val)
        assert str(c) == _pad("HISTORY " + val)
        val = "A commentary card has no comment."
        c = fits.Card("COMMENT", val, "comment")
        assert str(c) == _pad("COMMENT " + val)

    def test_commentary_card_created_by_fromstring(self):
        # commentary card created by fromstring()
        c = fits.Card.fromstring(
            "COMMENT card has no comments. "
            "/ text after slash is still part of the value."
        )
        assert (
            c.value == "card has no comments. "
            "/ text after slash is still part of the value."
        )
        assert c.comment == ""

    def test_commentary_card_will_not_parse_numerical_value(self):
        # commentary card will not parse the numerical value
        c = fits.Card.fromstring("HISTORY  (1, 2)")
        assert str(c) == _pad("HISTORY  (1, 2)")

    def test_equal_sign_after_column8(self):
        # equal sign after column 8 of a commentary card will be part of the
        # string value
        c = fits.Card.fromstring("HISTORY =   (1, 2)")
        assert str(c) == _pad("HISTORY =   (1, 2)")

    def test_blank_keyword(self):
        c = fits.Card("", "       / EXPOSURE INFORMATION")
        assert str(c) == _pad("               / EXPOSURE INFORMATION")
        c = fits.Card.fromstring(str(c))
        assert c.keyword == ""
        assert c.value == "       / EXPOSURE INFORMATION"

    def test_specify_undefined_value(self):
        # this is how to specify an undefined value
        c = fits.Card("undef", fits.card.UNDEFINED)
        assert str(c) == _pad("UNDEF   =")

    def test_complex_number_using_string_input(self):
        # complex number using string input
        c = fits.Card.fromstring("ABC     = (8, 9)")
        assert str(c) == _pad("ABC     = (8, 9)")

    def test_fixable_non_standard_fits_card(self, capsys):
        # fixable non-standard FITS card will keep the original format
        c = fits.Card.fromstring("abc     = +  2.1   e + 12")
        assert c.value == 2100000000000.0

        with pytest.warns(fits.verify.VerifyWarning) as w:
            assert str(c) == _pad("ABC     =             +2.1E+12")
        assert "Verification reported errors" in str(w[0].message)

    def test_fixable_non_fsc(self):
        # fixable non-FSC: if the card is not parsable, it's value will be
        # assumed
        # to be a string and everything after the first slash will be comment
        c = fits.Card.fromstring(
            "no_quote=  this card's value has no quotes / let's also try the comment"
        )

        with pytest.warns(fits.verify.VerifyWarning) as w:
            assert (
                str(c) == "NO_QUOTE= 'this card''s value has no quotes' "
                "/ let's also try the comment       "
            )
        assert "Verification reported errors" in str(w[0].message)

    def test_undefined_value_using_string_input(self):
        # undefined value using string input
        c = fits.Card.fromstring("ABC     =    ")
        assert str(c) == _pad("ABC     =")

    def test_undefined_keys_values(self):
        header = fits.Header()
        header["FOO"] = "BAR"
        header["UNDEF"] = None

        assert list(header.values()) == ["BAR", None]
        assert list(header.items()) == [("FOO", "BAR"), ("UNDEF", None)]

    def test_mislocated_equal_sign(self, capsys):
        # test mislocated "=" sign
        c = fits.Card.fromstring("XYZ= 100")
        assert c.keyword == "XYZ"
        assert c.value == 100

        with pytest.warns(fits.verify.VerifyWarning) as w:
            assert str(c) == _pad("XYZ     =                  100")
        assert "Verification reported errors" in str(w[0].message)

    def test_equal_only_up_to_column_10(self, capsys):
        # the test of "=" location is only up to column 10

        # This test used to check if Astropy rewrote this card to a new format,
        # something like "HISTO   = '=   (1, 2)".  But since ticket #109 if the
        # format is completely wrong we don't make any assumptions and the card
        # should be left alone
        c = fits.Card.fromstring("HISTO       =   (1, 2)")
        with pytest.warns(AstropyUserWarning, match=r"header keyword is invalid"):
            assert str(c) == _pad("HISTO       =   (1, 2)")

        # Likewise this card should just be left in its original form and
        # we shouldn't guess how to parse it or rewrite it.
        c = fits.Card.fromstring("   HISTORY          (1, 2)")
        with pytest.warns(AstropyUserWarning, match=r"header keyword is invalid"):
            assert str(c) == _pad("   HISTORY          (1, 2)")

    def test_verify_invalid_equal_sign(self):
        # verification
        c = fits.Card.fromstring("ABC= a6")
        with pytest.warns(AstropyUserWarning) as w:
            c.verify()
        err_text1 = "Card 'ABC' is not FITS standard (equal sign not at column 8)"
        err_text2 = "Card 'ABC' is not FITS standard (invalid value string: 'a6'"
        assert len(w) == 4
        assert err_text1 in str(w[1].message)
        assert err_text2 in str(w[2].message)

    def test_fix_invalid_equal_sign(self):
        fix_text = "Fixed 'ABC' card to meet the FITS standard."
        c = fits.Card.fromstring("ABC= a6")

        with pytest.warns(fits.verify.VerifyWarning) as w:
            c.verify("fix")
        assert len(w) == 4
        assert fix_text in str(w[2].message)
        assert str(c) == _pad("ABC     = 'a6      '")

    def test_long_string_value(self):
        # test long string value
        c = fits.Card("abc", "long string value " * 10, "long comment " * 10)
        assert (
            str(c)
            == "ABC     = 'long string value long string value long string value long string &' "
            "CONTINUE  'value long string value long string value long string value long &'  "
            "CONTINUE  'string value long string value long string value &'                  "
            "CONTINUE  '&' / long comment long comment long comment long comment long        "
            "CONTINUE  '&' / comment long comment long comment long comment long comment     "
            "CONTINUE  '' / long comment                                                     "
        )

    def test_long_string_value_with_multiple_long_words(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/11298
        """
        c = fits.Card(
            "WHATEVER",
            "SuperCalibrationParameters_XXXX_YYYY_ZZZZZ_KK_01_02_"
            "03)-AAABBBCCC.n.h5 SuperNavigationParameters_XXXX_YYYY"
            "_ZZZZZ_KK_01_02_03)-AAABBBCCC.n.xml",
        )
        assert (
            str(c)
            == "WHATEVER= 'SuperCalibrationParameters_XXXX_YYYY_ZZZZZ_KK_01_02_03)-AAABBBCCC.n&'"
            "CONTINUE  '.h5 &'                                                               "
            "CONTINUE  'SuperNavigationParameters_XXXX_YYYY_ZZZZZ_KK_01_02_03)-AAABBBCCC.n.&'"
            "CONTINUE  'xml'                                                                 "
        )

    def test_long_unicode_string(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/1

        So long as a unicode string can be converted to ASCII it should have no
        different behavior in this regard from a byte string.
        """

        h1 = fits.Header()
        h1["TEST"] = "abcdefg" * 30

        h2 = fits.Header()
        h2["TEST"] = "abcdefg" * 30

        assert str(h1) == str(h2)

    def test_long_string_repr(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/193

        Ensure that the __repr__() for cards represented with CONTINUE cards is
        split across multiple lines (broken at each *physical* card).
        """

        header = fits.Header()
        header["TEST1"] = ("Regular value", "Regular comment")
        header["TEST2"] = ("long string value " * 10, "long comment " * 10)
        header["TEST3"] = ("Regular value", "Regular comment")

        assert repr(header).splitlines() == [
            str(fits.Card("TEST1", "Regular value", "Regular comment")),
            "TEST2   = 'long string value long string value long string value long string &' ",
            "CONTINUE  'value long string value long string value long string value long &'  ",
            "CONTINUE  'string value long string value long string value &'                  ",
            "CONTINUE  '&' / long comment long comment long comment long comment long        ",
            "CONTINUE  '&' / comment long comment long comment long comment long comment     ",
            "CONTINUE  '' / long comment                                                     ",
            str(fits.Card("TEST3", "Regular value", "Regular comment")),
        ]

    def test_blank_keyword_long_value(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/194

        Test that a blank keyword ('') can be assigned a too-long value that is
        continued across multiple cards with blank keywords, just like COMMENT
        and HISTORY cards.
        """

        value = "long string value " * 10
        header = fits.Header()
        header[""] = value

        assert len(header) == 3
        assert " ".join(header[""]) == value.rstrip()

        # Ensure that this works like other commentary keywords
        header["COMMENT"] = value
        header["HISTORY"] = value
        assert header["COMMENT"] == header["HISTORY"]
        assert header["COMMENT"] == header[""]

    def test_long_string_from_file(self):
        c = fits.Card("abc", "long string value " * 10, "long comment " * 10)
        hdu = fits.PrimaryHDU()
        hdu.header.append(c)
        hdu.writeto(self.temp("test_new.fits"))

        hdul = fits.open(self.temp("test_new.fits"))
        c = hdul[0].header.cards["abc"]
        hdul.close()
        assert (
            str(c)
            == "ABC     = 'long string value long string value long string value long string &' "
            "CONTINUE  'value long string value long string value long string value long &'  "
            "CONTINUE  'string value long string value long string value &'                  "
            "CONTINUE  '&' / long comment long comment long comment long comment long        "
            "CONTINUE  '&' / comment long comment long comment long comment long comment     "
            "CONTINUE  '' / long comment                                                     "
        )

    def test_word_in_long_string_too_long(self):
        # if a word in a long string is too long, it will be cut in the middle
        c = fits.Card("abc", "longstringvalue" * 10, "longcomment" * 10)
        assert (
            str(c)
            == "ABC     = 'longstringvaluelongstringvaluelongstringvaluelongstringvaluelongstr&'"
            "CONTINUE  'ingvaluelongstringvaluelongstringvaluelongstringvaluelongstringvalu&'"
            "CONTINUE  'elongstringvalue&'                                                   "
            "CONTINUE  '&' / longcommentlongcommentlongcommentlongcommentlongcommentlongcomme"
            "CONTINUE  '' / ntlongcommentlongcommentlongcommentlongcomment                   "
        )

    def test_long_string_value_via_fromstring(self, capsys):
        # long string value via fromstring() method
        c = fits.Card.fromstring(
            _pad("abc     = 'longstring''s testing  &  ' / comments in line 1")
            + _pad(
                "continue  'continue with long string but without the "
                "ampersand at the end' /"
            )
            + _pad(
                "continue  'continue must have string value (with quotes)' "
                "/ comments with ''. "
            )
        )

        with pytest.warns(fits.verify.VerifyWarning) as w:
            assert (
                str(c)
                == "ABC     = 'longstring''s testing  continue with long string but without the &'  "
                "CONTINUE  'ampersand at the endcontinue must have string value (with quotes)&'  "
                "CONTINUE  '' / comments in line 1 comments with ''.                             "
            )
        assert "Verification reported errors" in str(w[0].message)

    def test_long_string_value_with_quotes(self):
        testval = "x" * 100 + "''"
        c = fits.Card("TEST", testval)
        c = fits.Card.fromstring(c.image)
        assert c.value == testval

        testval = "x" * 100 + "''xxx"
        c = fits.Card("TEST", testval)
        c = fits.Card.fromstring(c.image)
        assert c.value == testval

        testval = "x" * 100 + "'' xxx"
        c = fits.Card("TEST", testval)
        c = fits.Card.fromstring(c.image)
        assert c.value == testval

    def test_continue_card_with_equals_in_value(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/117
        """

        c = fits.Card.fromstring(
            _pad(
                "EXPR    = '/grp/hst/cdbs//grid/pickles/dat_uvk/pickles_uk_10.fits * &'"
            )
            + _pad("CONTINUE  '5.87359e-12 * MWAvg(Av=0.12)&'")
            + _pad("CONTINUE  '&' / pysyn expression")
        )

        assert c.keyword == "EXPR"
        assert (
            c.value == "/grp/hst/cdbs//grid/pickles/dat_uvk/pickles_uk_10.fits "
            "* 5.87359e-12 * MWAvg(Av=0.12)"
        )
        assert c.comment == "pysyn expression"

    def test_final_continue_card_lacks_ampersand(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/3282
        """

        h = fits.Header()
        h["SVALUE"] = "A" * 69
        assert repr(h).splitlines()[-1] == _pad("CONTINUE  'AA'")

    def test_final_continue_card_ampersand_removal_on_long_comments(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/3282
        """

        c = fits.Card("TEST", "long value" * 10, "long comment &" * 10)
        assert (
            str(c)
            == "TEST    = 'long valuelong valuelong valuelong valuelong valuelong valuelong &'  "
            "CONTINUE  'valuelong valuelong valuelong value&'                                "
            "CONTINUE  '&' / long comment &long comment &long comment &long comment &long    "
            "CONTINUE  '&' / comment &long comment &long comment &long comment &long comment "
            "CONTINUE  '' / &long comment &                                                  "
        )

    def test_hierarch_card_creation(self):
        # Test automatic upgrade to hierarch card
        with pytest.warns(
            AstropyUserWarning, match="HIERARCH card will be created"
        ) as w:
            c = fits.Card(
                "ESO INS SLIT2 Y1FRML",
                "ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)",
            )
        assert len(w) == 1
        assert (
            str(c) == "HIERARCH ESO INS SLIT2 Y1FRML= "
            "'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)'"
        )

        # Test manual creation of hierarch card
        c = fits.Card("hierarch abcdefghi", 10)
        assert str(c) == _pad("HIERARCH abcdefghi = 10")
        c = fits.Card(
            "HIERARCH ESO INS SLIT2 Y1FRML",
            "ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)",
        )
        assert (
            str(c) == "HIERARCH ESO INS SLIT2 Y1FRML= "
            "'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)'"
        )

    def test_hierarch_with_abbrev_value_indicator(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/5
        """

        c = fits.Card.fromstring("HIERARCH key.META_4='calFileVersion'")
        assert c.keyword == "key.META_4"
        assert c.value == "calFileVersion"
        assert c.comment == ""

    def test_hierarch_not_warn(self):
        """Check that compressed image headers do not issue HIERARCH warnings."""
        filename = fits.util.get_testdata_filepath("compressed_image.fits")
        with fits.open(filename) as hdul:
            header = hdul[1].header
        with warnings.catch_warnings(record=True) as warning_list:
            header["HIERARCH LONG KEYWORD"] = 42
        assert len(warning_list) == 0
        assert header["LONG KEYWORD"] == 42
        assert header["HIERARCH LONG KEYWORD"] == 42

        # Check that it still warns if we do not use HIERARCH
        with pytest.warns(
            fits.verify.VerifyWarning, match=r"greater than 8 characters"
        ):
            header["LONG KEYWORD2"] = 1
        assert header["LONG KEYWORD2"] == 1

    def test_hierarch_keyword_whitespace(self):
        """
        Regression test for
        https://github.com/spacetelescope/PyFITS/issues/6

        Make sure any leading or trailing whitespace around HIERARCH
        keywords is stripped from the actual keyword value.
        """

        c = fits.Card.fromstring("HIERARCH  key.META_4    = 'calFileVersion'")
        assert c.keyword == "key.META_4"
        assert c.value == "calFileVersion"
        assert c.comment == ""

        # Test also with creation via the Card constructor
        c = fits.Card("HIERARCH  key.META_4", "calFileVersion")
        assert c.keyword == "key.META_4"
        assert c.value == "calFileVersion"
        assert c.comment == ""

    def test_verify_mixed_case_hierarch(self):
        """Regression test for
        https://github.com/spacetelescope/PyFITS/issues/7

        Assures that HIERARCH keywords with lower-case characters and other
        normally invalid keyword characters are not considered invalid.
        """

        c = fits.Card("HIERARCH WeirdCard.~!@#_^$%&", "The value", "a comment")
        # This should not raise any exceptions
        c.verify("exception")
        assert c.keyword == "WeirdCard.~!@#_^$%&"
        assert c.value == "The value"
        assert c.comment == "a comment"

        # Test also the specific case from the original bug report
        header = fits.Header(
            [
                ("simple", True),
                ("BITPIX", 8),
                ("NAXIS", 0),
                ("EXTEND", True, "May contain datasets"),
                ("HIERARCH key.META_0", "detRow"),
            ]
        )
        hdu = fits.PrimaryHDU(header=header)
        hdu.writeto(self.temp("test.fits"))
        with fits.open(self.temp("test.fits")) as hdul:
            header2 = hdul[0].header
            assert str(header.cards[header.index("key.META_0")]) == str(
                header2.cards[header2.index("key.META_0")]
            )

    def test_missing_keyword(self):
        """Test that accessing a non-existent keyword raises a KeyError."""

        header = fits.Header()
        # De-referencing header through the inline function should behave
        # identically to accessing it in the pytest.raises context below.
        pytest.raises(KeyError, lambda k: header[k], "NAXIS")
        # Test exception with message
        with pytest.raises(KeyError, match=r"Keyword 'NAXIS' not found."):
            header["NAXIS"]

    def test_hierarch_card_lookup(self):
        header = fits.Header()
        header["hierarch abcdefghi"] = 10
        assert "abcdefghi" in header
        assert header["abcdefghi"] == 10
        # This used to be assert_false, but per ticket
        # https://aeon.stsci.edu/ssb/trac/pyfits/ticket/155 hierarch keywords
        # should be treated case-insensitively when performing lookups
        assert "ABCDEFGHI" in header

    def test_hierarch_card_delete(self):
        header = fits.Header()
        header["hierarch abcdefghi"] = 10
        del header["hierarch abcdefghi"]

    def test_hierarch_card_insert_delete(self):
        header = fits.Header()
        with pytest.warns(
            fits.verify.VerifyWarning, match=r"greater than 8 characters"
        ):
            header["abcdefghi"] = 10
        header["abcdefgh"] = 10
        header["abcdefg"] = 10
        with pytest.warns(
            fits.verify.VerifyWarning, match=r"greater than 8 characters"
        ):
            header.insert(2, ("abcdefghij", 10))
        del header["abcdefghij"]
        with pytest.warns(
            fits.verify.VerifyWarning, match=r"greater than 8 characters"
        ):
            header.insert(2, ("abcdefghij", 10))
        del header[2]
        assert list(header.keys())[2] == "abcdefg".upper()

    def test_hierarch_create_and_update(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/158

        Tests several additional use cases for working with HIERARCH cards.
        """

        msg = "a HIERARCH card will be created"

        header = fits.Header()
        with pytest.warns(VerifyWarning) as w:
            header.update({"HIERARCH BLAH BLAH": "TESTA"})
            assert len(w) == 0
            assert "BLAH BLAH" in header
            assert header["BLAH BLAH"] == "TESTA"

            header.update({"HIERARCH BLAH BLAH": "TESTB"})
            assert len(w) == 0
            assert header["BLAH BLAH"], "TESTB"

            # Update without explicitly stating 'HIERARCH':
            header.update({"BLAH BLAH": "TESTC"})
            assert len(w) == 1
            assert len(header) == 1
            assert header["BLAH BLAH"], "TESTC"

            # Test case-insensitivity
            header.update({"HIERARCH blah blah": "TESTD"})
            assert len(w) == 1
            assert len(header) == 1
            assert header["blah blah"], "TESTD"

            header.update({"blah blah": "TESTE"})
            assert len(w) == 2
            assert len(header) == 1
            assert header["blah blah"], "TESTE"

            # Create a HIERARCH card > 8 characters without explicitly stating
            # 'HIERARCH'
            header.update({"BLAH BLAH BLAH": "TESTA"})
            assert len(w) == 3
            assert msg in str(w[0].message)

            header.update({"HIERARCH BLAH BLAH BLAH": "TESTB"})
            assert len(w) == 3
            assert header["BLAH BLAH BLAH"], "TESTB"

            # Update without explicitly stating 'HIERARCH':
            header.update({"BLAH BLAH BLAH": "TESTC"})
            assert len(w) == 4
            assert header["BLAH BLAH BLAH"], "TESTC"

            # Test case-insensitivity
            header.update({"HIERARCH blah blah blah": "TESTD"})
            assert len(w) == 4
            assert header["blah blah blah"], "TESTD"

            header.update({"blah blah blah": "TESTE"})
            assert len(w) == 5
            assert header["blah blah blah"], "TESTE"

    def test_short_hierarch_create_and_update(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/158

        Tests several additional use cases for working with HIERARCH cards,
        specifically where the keyword is fewer than 8 characters, but contains
        invalid characters such that it can only be created as a HIERARCH card.
        """

        msg = "a HIERARCH card will be created"

        header = fits.Header()
        with pytest.warns(VerifyWarning) as w:
            header.update({"HIERARCH BLA BLA": "TESTA"})
            assert len(w) == 0
            assert "BLA BLA" in header
            assert header["BLA BLA"] == "TESTA"

            header.update({"HIERARCH BLA BLA": "TESTB"})
            assert len(w) == 0
            assert header["BLA BLA"], "TESTB"

            # Update without explicitly stating 'HIERARCH':
            header.update({"BLA BLA": "TESTC"})
            assert len(w) == 1
            assert header["BLA BLA"], "TESTC"

            # Test case-insensitivity
            header.update({"HIERARCH bla bla": "TESTD"})
            assert len(w) == 1
            assert len(header) == 1
            assert header["bla bla"], "TESTD"

            header.update({"bla bla": "TESTE"})
            assert len(w) == 2
            assert len(header) == 1
            assert header["bla bla"], "TESTE"

        header = fits.Header()
        with pytest.warns(VerifyWarning) as w:
            # Create a HIERARCH card containing invalid characters without
            # explicitly stating 'HIERARCH'
            header.update({"BLA BLA": "TESTA"})
            print([x.category for x in w])
            assert len(w) == 1
            assert msg in str(w[0].message)

            header.update({"HIERARCH BLA BLA": "TESTB"})
            assert len(w) == 1
            assert header["BLA BLA"], "TESTB"

            # Update without explicitly stating 'HIERARCH':
            header.update({"BLA BLA": "TESTC"})
            assert len(w) == 2
            assert header["BLA BLA"], "TESTC"

            # Test case-insensitivity
            header.update({"HIERARCH bla bla": "TESTD"})
            assert len(w) == 2
            assert len(header) == 1
            assert header["bla bla"], "TESTD"

            header.update({"bla bla": "TESTE"})
            assert len(w) == 3
            assert len(header) == 1
            assert header["bla bla"], "TESTE"

    def test_header_setitem_invalid(self):
        header = fits.Header()

        def test():
            header["FOO"] = ("bar", "baz", "qux")

        pytest.raises(ValueError, test)

    def test_header_setitem_1tuple(self):
        header = fits.Header()
        header["FOO"] = ("BAR",)
        header["FOO2"] = (None,)
        assert header["FOO"] == "BAR"
        assert header["FOO2"] is None
        assert header[0] == "BAR"
        assert header.comments[0] == ""
        assert header.comments["FOO"] == ""

    def test_header_setitem_2tuple(self):
        header = fits.Header()
        header["FOO"] = ("BAR", "BAZ")
        header["FOO2"] = (None, None)
        assert header["FOO"] == "BAR"
        assert header["FOO2"] is None
        assert header[0] == "BAR"
        assert header.comments[0] == "BAZ"
        assert header.comments["FOO"] == "BAZ"
        assert header.comments["FOO2"] == ""

    def test_header_set_value_to_none(self):
        """
        Setting the value of a card to None should simply give that card an
        undefined value.  Undefined value should map to None.
        """

        header = fits.Header()
        header["FOO"] = "BAR"
        assert header["FOO"] == "BAR"
        header["FOO"] = None
        assert header["FOO"] is None

        # Create a header that contains an undefined value and a defined
        # value.
        hstr = "UNDEF   = \nDEFINED = 42"
        header = fits.Header.fromstring(hstr, sep="\n")

        # Explicitly add a card with an UNDEFINED value
        c = fits.Card("UNDEF2", fits.card.UNDEFINED)
        header.extend([c])

        # And now assign an undefined value to the header through setitem
        header["UNDEF3"] = fits.card.UNDEFINED

        # Tuple assignment
        header.append(("UNDEF5", None, "Undefined value"), end=True)
        header.append("UNDEF6")

        assert header["DEFINED"] == 42
        assert header["UNDEF"] is None
        assert header["UNDEF2"] is None
        assert header["UNDEF3"] is None
        assert header["UNDEF5"] is None
        assert header["UNDEF6"] is None

        # Assign an undefined value to a new card
        header["UNDEF4"] = None

        # Overwrite an existing value with None
        header["DEFINED"] = None

        # All headers now should be undefined
        for c in header.cards:
            assert c.value == fits.card.UNDEFINED

    def test_set_comment_only(self):
        header = fits.Header([("A", "B", "C")])
        header.set("A", comment="D")
        assert header["A"] == "B"
        assert header.comments["A"] == "D"

    def test_header_iter(self):
        header = fits.Header([("A", "B"), ("C", "D")])
        assert list(header) == ["A", "C"]

    def test_header_slice(self):
        header = fits.Header([("A", "B"), ("C", "D"), ("E", "F")])
        newheader = header[1:]
        assert len(newheader) == 2
        assert "A" not in newheader
        assert "C" in newheader
        assert "E" in newheader

        newheader = header[::-1]
        assert len(newheader) == 3
        assert newheader[0] == "F"
        assert newheader[1] == "D"
        assert newheader[2] == "B"

        newheader = header[::2]
        assert len(newheader) == 2
        assert "A" in newheader
        assert "C" not in newheader
        assert "E" in newheader

    def test_header_slice_assignment(self):
        """
        Assigning to a slice should just assign new values to the cards
        included in the slice.
        """

        header = fits.Header([("A", "B"), ("C", "D"), ("E", "F")])

        # Test assigning slice to the same value; this works similarly to numpy
        # arrays
        header[1:] = 1
        assert header[1] == 1
        assert header[2] == 1

        # Though strings are iterable they should be treated as a scalar value
        header[1:] = "GH"
        assert header[1] == "GH"
        assert header[2] == "GH"

        # Now assign via an iterable
        header[1:] = ["H", "I"]
        assert header[1] == "H"
        assert header[2] == "I"

    def test_header_slice_delete(self):
        """Test deleting a slice of cards from the header."""

        header = fits.Header([("A", "B"), ("C", "D"), ("E", "F")])
        del header[1:]
        assert len(header) == 1
        assert header[0] == "B"
        del header[:]
        assert len(header) == 0

    def test_wildcard_slice(self):
        """Test selecting a subsection of a header via wildcard matching."""

        header = fits.Header([("ABC", 0), ("DEF", 1), ("ABD", 2)])
        newheader = header["AB*"]
        assert len(newheader) == 2
        assert newheader[0] == 0
        assert newheader[1] == 2

    def test_wildcard_with_hyphen(self):
        """
        Regression test for issue where wildcards did not work on keywords
        containing hyphens.
        """

        header = fits.Header([("DATE", 1), ("DATE-OBS", 2), ("DATE-FOO", 3)])
        assert len(header["DATE*"]) == 3
        assert len(header["DATE?*"]) == 2
        assert len(header["DATE-*"]) == 2

    def test_wildcard_slice_assignment(self):
        """Test assigning to a header slice selected via wildcard matching."""

        header = fits.Header([("ABC", 0), ("DEF", 1), ("ABD", 2)])

        # Test assigning slice to the same value; this works similarly to numpy
        # arrays
        header["AB*"] = 1
        assert header[0] == 1
        assert header[2] == 1

        # Though strings are iterable they should be treated as a scalar value
        header["AB*"] = "GH"
        assert header[0] == "GH"
        assert header[2] == "GH"

        # Now assign via an iterable
        header["AB*"] = ["H", "I"]
        assert header[0] == "H"
        assert header[2] == "I"

    def test_wildcard_slice_deletion(self):
        """Test deleting cards from a header that match a wildcard pattern."""

        header = fits.Header([("ABC", 0), ("DEF", 1), ("ABD", 2)])
        del header["AB*"]
        assert len(header) == 1
        assert header[0] == 1

    def test_header_history(self):
        header = fits.Header(
            [
                ("ABC", 0),
                ("HISTORY", 1),
                ("HISTORY", 2),
                ("DEF", 3),
                ("HISTORY", 4),
                ("HISTORY", 5),
            ]
        )
        assert header["HISTORY"] == [1, 2, 4, 5]

    def test_header_clear(self):
        header = fits.Header([("A", "B"), ("C", "D")])
        header.clear()
        assert "A" not in header
        assert "C" not in header
        assert len(header) == 0

    @pytest.mark.parametrize("fitsext", [fits.ImageHDU(), fits.CompImageHDU()])
    def test_header_clear_write(self, fitsext):
        hdulist = fits.HDUList([fits.PrimaryHDU(), fitsext])
        hdulist[1].header["FOO"] = "BAR"
        hdulist[1].header.clear()
        with pytest.raises(VerifyError) as err:
            hdulist.writeto(self.temp("temp.fits"), overwrite=True)
        err_msg = "'XTENSION' card does not exist."
        assert err_msg in str(err.value)

    def test_header_fromkeys(self):
        header = fits.Header.fromkeys(["A", "B"])
        assert "A" in header
        assert header["A"] is None
        assert header.comments["A"] == ""
        assert "B" in header
        assert header["B"] is None
        assert header.comments["B"] == ""

    def test_header_fromkeys_with_value(self):
        header = fits.Header.fromkeys(["A", "B"], "C")
        assert "A" in header
        assert header["A"] == "C"
        assert header.comments["A"] == ""
        assert "B" in header
        assert header["B"] == "C"
        assert header.comments["B"] == ""

    def test_header_fromkeys_with_value_and_comment(self):
        header = fits.Header.fromkeys(["A"], ("B", "C"))
        assert "A" in header
        assert header["A"] == "B"
        assert header.comments["A"] == "C"

    def test_header_fromkeys_with_duplicates(self):
        header = fits.Header.fromkeys(["A", "B", "A"], "C")
        assert "A" in header
        assert ("A", 0) in header
        assert ("A", 1) in header
        assert ("A", 2) not in header
        assert header[0] == "C"
        assert header["A"] == "C"
        assert header[("A", 0)] == "C"
        assert header[2] == "C"
        assert header[("A", 1)] == "C"

    def test_header_items(self):
        header = fits.Header([("A", "B"), ("C", "D")])
        assert list(header.items()) == [("A", "B"), ("C", "D")]

    def test_header_iterkeys(self):
        header = fits.Header([("A", "B"), ("C", "D")])
        for a, b in zip(header.keys(), header):
            assert a == b

    def test_header_itervalues(self):
        header = fits.Header([("A", "B"), ("C", "D")])
        for a, b in zip(header.values(), ["B", "D"]):
            assert a == b

    def test_header_keys(self):
        with fits.open(self.data("arange.fits")) as hdul:
            assert list(hdul[0].header) == [
                "SIMPLE",
                "BITPIX",
                "NAXIS",
                "NAXIS1",
                "NAXIS2",
                "NAXIS3",
                "EXTEND",
            ]

    def test_header_list_like_pop(self):
        header = fits.Header([("A", "B"), ("C", "D"), ("E", "F"), ("G", "H")])

        last = header.pop()
        assert last == "H"
        assert len(header) == 3
        assert list(header) == ["A", "C", "E"]

        mid = header.pop(1)
        assert mid == "D"
        assert len(header) == 2
        assert list(header) == ["A", "E"]

        first = header.pop(0)
        assert first == "B"
        assert len(header) == 1
        assert list(header) == ["E"]

        pytest.raises(IndexError, header.pop, 42)

    def test_header_dict_like_pop(self):
        header = fits.Header([("A", "B"), ("C", "D"), ("E", "F"), ("G", "H")])
        pytest.raises(TypeError, header.pop, "A", "B", "C")

        last = header.pop("G")
        assert last == "H"
        assert len(header) == 3
        assert list(header) == ["A", "C", "E"]

        mid = header.pop("C")
        assert mid == "D"
        assert len(header) == 2
        assert list(header) == ["A", "E"]

        first = header.pop("A")
        assert first == "B"
        assert len(header) == 1
        assert list(header) == ["E"]

        default = header.pop("X", "Y")
        assert default == "Y"
        assert len(header) == 1

        pytest.raises(KeyError, header.pop, "X")

    def test_popitem(self):
        header = fits.Header([("A", "B"), ("C", "D"), ("E", "F")])
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
        header = fits.Header([("A", "B"), ("C", "D"), ("E", "F")])
        assert header.setdefault("A") == "B"
        assert header.setdefault("C") == "D"
        assert header.setdefault("E") == "F"
        assert len(header) == 3
        assert header.setdefault("G", "H") == "H"
        assert len(header) == 4
        assert "G" in header
        assert header.setdefault("G", "H") == "H"
        assert len(header) == 4

    def test_update_from_dict(self):
        """
        Test adding new cards and updating existing cards from a dict using
        Header.update()
        """

        header = fits.Header([("A", "B"), ("C", "D")])
        header.update({"A": "E", "F": "G"})
        assert header["A"] == "E"
        assert header[0] == "E"
        assert "F" in header
        assert header["F"] == "G"
        assert header[-1] == "G"

        # Same as above but this time pass the update dict as keyword arguments
        header = fits.Header([("A", "B"), ("C", "D")])
        header.update(A="E", F="G")
        assert header["A"] == "E"
        assert header[0] == "E"
        assert "F" in header
        assert header["F"] == "G"
        assert header[-1] == "G"

    def test_update_from_iterable(self):
        """
        Test adding new cards and updating existing cards from an iterable of
        cards and card tuples.
        """

        header = fits.Header([("A", "B"), ("C", "D")])
        header.update([("A", "E"), fits.Card("F", "G")])
        assert header["A"] == "E"
        assert header[0] == "E"
        assert "F" in header
        assert header["F"] == "G"
        assert header[-1] == "G"

    def test_header_extend(self):
        """
        Test extending a header both with and without stripping cards from the
        extension header.
        """

        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu2.header["MYKEY"] = ("some val", "some comment")
        hdu.header += hdu2.header
        assert len(hdu.header) == 5
        assert hdu.header[-1] == "some val"

        # Same thing, but using + instead of +=
        hdu = fits.PrimaryHDU()
        hdu.header = hdu.header + hdu2.header
        assert len(hdu.header) == 5
        assert hdu.header[-1] == "some val"

        # Directly append the other header in full--not usually a desirable
        # operation when the header is coming from another HDU
        hdu.header.extend(hdu2.header, strip=False)
        assert len(hdu.header) == 11
        assert list(hdu.header)[5] == "XTENSION"
        assert hdu.header[-1] == "some val"
        assert ("MYKEY", 1) in hdu.header

    def test_header_extend_unique(self):
        """
        Test extending the header with and without unique=True.
        """
        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu.header["MYKEY"] = ("some val", "some comment")
        hdu2.header["MYKEY"] = ("some other val", "some other comment")
        hdu.header.extend(hdu2.header)
        assert len(hdu.header) == 6
        assert hdu.header[-2] == "some val"
        assert hdu.header[-1] == "some other val"

        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu.header["MYKEY"] = ("some val", "some comment")
        hdu2.header["MYKEY"] = ("some other val", "some other comment")
        hdu.header.extend(hdu2.header, unique=True)
        assert len(hdu.header) == 5
        assert hdu.header[-1] == "some val"

    def test_header_extend_unique_commentary(self):
        """
        Test extending header with and without unique=True and commentary
        cards in the header being added. Issue astropy/astropy#3967
        """
        for commentary_card in ["", "COMMENT", "HISTORY"]:
            for is_unique in [True, False]:
                hdu = fits.PrimaryHDU()
                # Make sure we are testing the case we want.
                assert commentary_card not in hdu.header
                hdu2 = fits.ImageHDU()
                hdu2.header[commentary_card] = "My text"
                hdu.header.extend(hdu2.header, unique=is_unique)
                assert len(hdu.header) == 5
                assert hdu.header[commentary_card][0] == "My text"

    def test_header_extend_update(self):
        """
        Test extending the header with and without update=True.
        """

        hdu = fits.PrimaryHDU()
        hdu2 = fits.ImageHDU()
        hdu.header["MYKEY"] = ("some val", "some comment")
        hdu.header["HISTORY"] = "history 1"
        hdu2.header["MYKEY"] = ("some other val", "some other comment")
        hdu2.header["HISTORY"] = "history 1"
        hdu2.header["HISTORY"] = "history 2"
        hdu.header.extend(hdu2.header)
        assert len(hdu.header) == 9
        assert ("MYKEY", 0) in hdu.header
        assert ("MYKEY", 1) in hdu.header
        assert hdu.header[("MYKEY", 1)] == "some other val"
        assert len(hdu.header["HISTORY"]) == 3
        assert hdu.header[-1] == "history 2"

        hdu = fits.PrimaryHDU()
        hdu.header["MYKEY"] = ("some val", "some comment")
        hdu.header["HISTORY"] = "history 1"
        hdu.header.extend(hdu2.header, update=True)
        assert len(hdu.header) == 7
        assert ("MYKEY", 0) in hdu.header
        assert ("MYKEY", 1) not in hdu.header
        assert hdu.header["MYKEY"] == "some other val"
        assert len(hdu.header["HISTORY"]) == 2
        assert hdu.header[-1] == "history 2"

    def test_header_extend_update_commentary(self):
        """
        Test extending header with and without unique=True and commentary
        cards in the header being added.

        Though not quite the same as astropy/astropy#3967, update=True hits
        the same if statement as that issue.
        """
        for commentary_card in ["", "COMMENT", "HISTORY"]:
            for is_update in [True, False]:
                hdu = fits.PrimaryHDU()
                # Make sure we are testing the case we want.
                assert commentary_card not in hdu.header
                hdu2 = fits.ImageHDU()
                hdu2.header[commentary_card] = "My text"
                hdu.header.extend(hdu2.header, update=is_update)
                assert len(hdu.header) == 5
                assert hdu.header[commentary_card][0] == "My text"

    def test_header_extend_exact(self):
        """
        Test that extending an empty header with the contents of an existing
        header can exactly duplicate that header, given strip=False and
        end=True.
        """

        header = fits.getheader(self.data("test0.fits"))
        header2 = fits.Header()
        header2.extend(header, strip=False, end=True)
        assert header == header2

    def test_header_count(self):
        header = fits.Header([("A", "B"), ("C", "D"), ("E", "F")])
        assert header.count("A") == 1
        assert header.count("C") == 1
        assert header.count("E") == 1
        header["HISTORY"] = "a"
        header["HISTORY"] = "b"
        assert header.count("HISTORY") == 2
        pytest.raises(KeyError, header.count, "G")

    def test_header_append_use_blanks(self):
        """
        Tests that blank cards can be appended, and that future appends will
        use blank cards when available (unless useblanks=False)
        """

        header = fits.Header([("A", "B"), ("C", "D")])

        # Append a couple blanks
        header.append()
        header.append()
        assert len(header) == 4
        assert header[-1] == ""
        assert header[-2] == ""

        # New card should fill the first blank by default
        header.append(("E", "F"))
        assert len(header) == 4
        assert header[-2] == "F"
        assert header[-1] == ""

        # This card should not use up a blank spot
        header.append(("G", "H"), useblanks=False)
        assert len(header) == 5
        assert header[-1] == ""
        assert header[-2] == "H"

    def test_header_append_keyword_only(self):
        """
        Test appending a new card with just the keyword, and no value or
        comment given.
        """

        header = fits.Header([("A", "B"), ("C", "D")])
        header.append("E")
        assert len(header) == 3
        assert list(header)[-1] == "E"
        assert header[-1] is None
        assert header.comments["E"] == ""

        # Try appending a blank--normally this can be accomplished with just
        # header.append(), but header.append('') should also work (and is maybe
        # a little more clear)
        header.append("")
        assert len(header) == 4

        assert list(header)[-1] == ""
        assert header[""] == ""
        assert header.comments[""] == ""

    def test_header_insert_use_blanks(self):
        header = fits.Header([("A", "B"), ("C", "D")])

        # Append a couple blanks
        header.append()
        header.append()

        # Insert a new card; should use up one of the blanks
        header.insert(1, ("E", "F"))
        assert len(header) == 4
        assert header[1] == "F"
        assert header[-1] == ""
        assert header[-2] == "D"

        # Insert a new card without using blanks
        header.insert(1, ("G", "H"), useblanks=False)
        assert len(header) == 5
        assert header[1] == "H"
        assert header[-1] == ""

    def test_header_insert_before_keyword(self):
        """
        Test that a keyword name or tuple can be used to insert new keywords.

        Also tests the ``after`` keyword argument.

        Regression test for https://github.com/spacetelescope/PyFITS/issues/12
        """

        header = fits.Header(
            [("NAXIS1", 10), ("COMMENT", "Comment 1"), ("COMMENT", "Comment 3")]
        )

        header.insert("NAXIS1", ("NAXIS", 2, "Number of axes"))
        assert list(header.keys())[0] == "NAXIS"
        assert header[0] == 2
        assert header.comments[0] == "Number of axes"

        header.insert("NAXIS1", ("NAXIS2", 20), after=True)
        assert list(header.keys())[1] == "NAXIS1"
        assert list(header.keys())[2] == "NAXIS2"
        assert header[2] == 20

        header.insert(("COMMENT", 1), ("COMMENT", "Comment 2"))
        assert header["COMMENT"] == ["Comment 1", "Comment 2", "Comment 3"]

        header.insert(("COMMENT", 2), ("COMMENT", "Comment 4"), after=True)
        assert header["COMMENT"] == ["Comment 1", "Comment 2", "Comment 3", "Comment 4"]

        header.insert(-1, ("TEST1", True))
        assert list(header.keys())[-2] == "TEST1"

        header.insert(-1, ("TEST2", True), after=True)
        assert list(header.keys())[-1] == "TEST2"
        assert list(header.keys())[-3] == "TEST1"

    def test_remove(self):
        header = fits.Header([("A", "B"), ("C", "D")])

        # When keyword is present in the header it should be removed.
        header.remove("C")
        assert len(header) == 1
        assert list(header) == ["A"]
        assert "C" not in header

        # When keyword is not present in the header and ignore_missing is
        # False, KeyError should be raised
        with pytest.raises(KeyError):
            header.remove("F")

        # When keyword is not present and ignore_missing is True, KeyError
        # will be ignored
        header.remove("F", ignore_missing=True)
        assert len(header) == 1

        # Test for removing all instances of a keyword
        header = fits.Header([("A", "B"), ("C", "D"), ("A", "F")])
        header.remove("A", remove_all=True)
        assert "A" not in header
        assert len(header) == 1
        assert list(header) == ["C"]
        assert header[0] == "D"

    def test_header_comments(self):
        header = fits.Header([("A", "B", "C"), ("DEF", "G", "H")])
        assert repr(header.comments) == "       A  C\n     DEF  H"

    def test_comment_slices_and_filters(self):
        header = fits.Header([("AB", "C", "D"), ("EF", "G", "H"), ("AI", "J", "K")])
        s = header.comments[1:]
        assert list(s) == ["H", "K"]
        s = header.comments[::-1]
        assert list(s) == ["K", "H", "D"]
        s = header.comments["A*"]
        assert list(s) == ["D", "K"]

    def test_comment_slice_filter_assign(self):
        header = fits.Header([("AB", "C", "D"), ("EF", "G", "H"), ("AI", "J", "K")])
        header.comments[1:] = "L"
        assert list(header.comments) == ["D", "L", "L"]
        assert header.cards[header.index("AB")].comment == "D"
        assert header.cards[header.index("EF")].comment == "L"
        assert header.cards[header.index("AI")].comment == "L"

        header.comments[::-1] = header.comments[:]
        assert list(header.comments) == ["L", "L", "D"]

        header.comments["A*"] = ["M", "N"]
        assert list(header.comments) == ["M", "L", "N"]

    def test_commentary_slicing(self):
        header = fits.Header()

        indices = list(range(5))

        for idx in indices:
            header["HISTORY"] = idx

        # Just a few sample slice types; this won't get all corner cases but if
        # these all work we should be in good shape
        assert header["HISTORY"][1:] == indices[1:]
        assert header["HISTORY"][:3] == indices[:3]
        assert header["HISTORY"][:6] == indices[:6]
        assert header["HISTORY"][:-2] == indices[:-2]
        assert header["HISTORY"][::-1] == indices[::-1]
        assert header["HISTORY"][1::-1] == indices[1::-1]
        assert header["HISTORY"][1:5:2] == indices[1:5:2]

        # Same tests, but copy the values first; as it turns out this is
        # different from just directly doing an __eq__ as in the first set of
        # assertions
        header.insert(0, ("A", "B", "C"))
        header.append(("D", "E", "F"), end=True)
        assert list(header["HISTORY"][1:]) == indices[1:]
        assert list(header["HISTORY"][:3]) == indices[:3]
        assert list(header["HISTORY"][:6]) == indices[:6]
        assert list(header["HISTORY"][:-2]) == indices[:-2]
        assert list(header["HISTORY"][::-1]) == indices[::-1]
        assert list(header["HISTORY"][1::-1]) == indices[1::-1]
        assert list(header["HISTORY"][1:5:2]) == indices[1:5:2]

    def test_update_commentary(self):
        header = fits.Header()
        header["FOO"] = "BAR"
        header["HISTORY"] = "ABC"
        header["FRED"] = "BARNEY"
        header["HISTORY"] = "DEF"
        header["HISTORY"] = "GHI"

        assert header["HISTORY"] == ["ABC", "DEF", "GHI"]

        # Single value update
        header["HISTORY"][0] = "FOO"
        assert header["HISTORY"] == ["FOO", "DEF", "GHI"]

        # Single value partial slice update
        header["HISTORY"][1:] = "BAR"
        assert header["HISTORY"] == ["FOO", "BAR", "BAR"]

        # Multi-value update
        header["HISTORY"][:] = ["BAZ", "QUX"]
        assert header["HISTORY"] == ["BAZ", "QUX", "BAR"]

    def test_commentary_comparison(self):
        """
        Regression test for an issue found in *writing* the regression test for
        https://github.com/astropy/astropy/issues/2363, where comparison of
        the list of values for a commentary keyword did not always compare
        correctly with other iterables.
        """

        header = fits.Header()
        header["HISTORY"] = "hello world"
        header["HISTORY"] = "hello world"
        header["COMMENT"] = "hello world"
        assert header["HISTORY"] != header["COMMENT"]
        header["COMMENT"] = "hello world"
        assert header["HISTORY"] == header["COMMENT"]

    def test_long_commentary_card(self):
        header = fits.Header()
        header["FOO"] = "BAR"
        header["BAZ"] = "QUX"
        longval = "ABC" * 30
        header["HISTORY"] = longval
        header["FRED"] = "BARNEY"
        header["HISTORY"] = longval

        assert len(header) == 7
        assert list(header)[2] == "FRED"
        assert str(header.cards[3]) == "HISTORY " + longval[:72]
        assert str(header.cards[4]).rstrip() == "HISTORY " + longval[72:]

        header.set("HISTORY", longval, after="FOO")
        assert len(header) == 9
        assert str(header.cards[1]) == "HISTORY " + longval[:72]
        assert str(header.cards[2]).rstrip() == "HISTORY " + longval[72:]

        header = fits.Header()
        header.update({"FOO": "BAR"})
        header.update({"BAZ": "QUX"})
        longval = "ABC" * 30
        header.add_history(longval)
        header.update({"FRED": "BARNEY"})
        header.add_history(longval)

        assert len(header.cards) == 7
        assert header.cards[2].keyword == "FRED"
        assert str(header.cards[3]) == "HISTORY " + longval[:72]
        assert str(header.cards[4]).rstrip() == "HISTORY " + longval[72:]

        header.add_history(longval, after="FOO")
        assert len(header.cards) == 9
        assert str(header.cards[1]) == "HISTORY " + longval[:72]
        assert str(header.cards[2]).rstrip() == "HISTORY " + longval[72:]

    def test_totxtfile(self, home_is_temp):
        header_filename = self.temp("header.txt")
        with fits.open(self.data("test0.fits")) as hdul:
            hdul[0].header.totextfile(header_filename)
            # Check the `overwrite` flag
            with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
                hdul[0].header.totextfile(header_filename, overwrite=False)
            hdul[0].header.totextfile(header_filename, overwrite=True)

        hdu = fits.ImageHDU()
        hdu.header.update({"MYKEY": "FOO"})
        hdu.header.extend(
            hdu.header.fromtextfile(header_filename), update=True, update_first=True
        )

        # Write the hdu out and read it back in again--it should be recognized
        # as a PrimaryHDU
        hdu.writeto(self.temp("test.fits"), output_verify="ignore")

        with fits.open(self.temp("test.fits")) as hdul:
            assert isinstance(hdul[0], fits.PrimaryHDU)

        hdu = fits.ImageHDU()
        hdu.header.update({"MYKEY": "FOO"})
        hdu.header.extend(
            hdu.header.fromtextfile(header_filename),
            update=True,
            update_first=True,
            strip=False,
        )
        assert "MYKEY" in hdu.header
        assert "EXTENSION" not in hdu.header
        assert "SIMPLE" in hdu.header

        hdu.writeto(self.temp("test.fits"), output_verify="ignore", overwrite=True)

        with fits.open(self.temp("test.fits")) as hdul2:
            assert len(hdul2) == 2
            assert "MYKEY" in hdul2[1].header

    def test_tofile(self, home_is_temp):
        """
        Repeat test_totxtfile, but with tofile()
        """
        header_filename = self.temp("header.fits")
        with fits.open(self.data("test0.fits")) as hdul:
            hdul[0].header.tofile(header_filename)
            # Check the `overwrite` flag
            with pytest.raises(OSError, match=_NOT_OVERWRITING_MSG_MATCH):
                hdul[0].header.tofile(header_filename, overwrite=False)
            hdul[0].header.tofile(header_filename, overwrite=True)

        hdu = fits.ImageHDU()
        hdu.header.update({"MYKEY": "FOO"})
        hdu.header.extend(
            hdu.header.fromfile(header_filename), update=True, update_first=True
        )

        # Write the hdu out and read it back in again--it should be recognized
        # as a PrimaryHDU
        hdu.writeto(self.temp("test.fits"), output_verify="ignore")

        with fits.open(self.temp("test.fits")) as hdul:
            assert isinstance(hdul[0], fits.PrimaryHDU)

        hdu = fits.ImageHDU()
        hdu.header.update({"MYKEY": "FOO"})
        hdu.header.extend(
            hdu.header.fromfile(header_filename),
            update=True,
            update_first=True,
            strip=False,
        )
        assert "MYKEY" in hdu.header
        assert "EXTENSION" not in hdu.header
        assert "SIMPLE" in hdu.header

        hdu.writeto(self.temp("test.fits"), output_verify="ignore", overwrite=True)

        with fits.open(self.temp("test.fits")) as hdul2:
            assert len(hdul2) == 2
            assert "MYKEY" in hdul2[1].header

    def test_fromfile(self):
        """Regression test for https://github.com/astropy/astropy/issues/8711"""
        filename = self.data("scale.fits")
        hdr = fits.Header.fromfile(filename)
        assert hdr["DATASET"] == "2MASS"

    def test_header_fromtextfile(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/122

        Manually write a text file containing some header cards ending with
        newlines and ensure that fromtextfile can read them back in.
        """

        header = fits.Header()
        header["A"] = ("B", "C")
        header["B"] = ("C", "D")
        header["C"] = ("D", "E")

        with open(self.temp("test.hdr"), "w") as f:
            f.write("\n".join(str(c).strip() for c in header.cards))

        header2 = fits.Header.fromtextfile(self.temp("test.hdr"))
        assert header == header2

    def test_header_fromtextfile_with_end_card(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/154

        Make sure that when a Header is read from a text file that the END card
        is ignored.
        """

        header = fits.Header([("A", "B", "C"), ("D", "E", "F")])

        # We don't use header.totextfile here because it writes each card with
        # trailing spaces to pad them out to 80 characters.  But this bug only
        # presents itself when each card ends immediately with a newline, and
        # no trailing spaces
        with open(self.temp("test.hdr"), "w") as f:
            f.write("\n".join(str(c).strip() for c in header.cards))
            f.write("\nEND")

        new_header = fits.Header.fromtextfile(self.temp("test.hdr"))

        assert "END" not in new_header
        assert header == new_header

    def test_append_end_card(self):
        """
        Regression test 2 for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/154

        Manually adding an END card to a header should simply result in a
        ValueError (as was the case in PyFITS 3.0 and earlier).
        """

        header = fits.Header([("A", "B", "C"), ("D", "E", "F")])

        def setitem(k, v):
            header[k] = v

        pytest.raises(ValueError, setitem, "END", "")
        pytest.raises(ValueError, header.append, "END")
        pytest.raises(ValueError, header.append, "END", end=True)
        pytest.raises(ValueError, header.insert, len(header), "END")
        pytest.raises(ValueError, header.set, "END")

    def test_invalid_end_cards(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/217

        This tests the case where the END card looks like a normal card like
        'END = ' and other similar oddities.  As long as a card starts with END
        and looks like it was intended to be the END card we allow it, but with
        a warning.
        """

        horig = fits.PrimaryHDU(data=np.arange(100)).header

        def invalid_header(end, pad):
            # Build up a goofy invalid header
            # Start from a seemingly normal header
            s = horig.tostring(sep="", endcard=False, padding=False)
            # append the bogus end card
            s += end
            # add additional padding if requested
            if pad:
                s += " " * _pad_length(len(s))

            # This will differ between Python versions
            if isinstance(s, bytes):
                return BytesIO(s)
            else:
                return StringIO(s)

        # Basic case motivated by the original issue; it's as if the END card
        # was appended by software that doesn't know to treat it specially, and
        # it is given an = after it
        s = invalid_header("END =", True)

        with pytest.warns(
            AstropyUserWarning, match="Unexpected bytes trailing END keyword: ' ='"
        ) as w:
            h = fits.Header.fromfile(s)
            assert h == horig
        assert len(w) == 1

        # A case similar to the last but with more spaces between END and the
        # =, as though the '= ' value indicator were placed like that of a
        # normal card
        s = invalid_header("END     = ", True)
        with pytest.warns(
            AstropyUserWarning, match="Unexpected bytes trailing END keyword: '     ='"
        ) as w:
            h = fits.Header.fromfile(s)
            assert h == horig
        assert len(w) == 1

        # END card with trailing gibberish
        s = invalid_header("END$%&%^*%*", True)
        with pytest.warns(
            AstropyUserWarning,
            match=r"Unexpected bytes trailing END keyword: '\$%&%\^\*%\*'",
        ) as w:
            h = fits.Header.fromfile(s)
            assert h == horig
        assert len(w) == 1

        # 'END' at the very end of a truncated file without padding; the way
        # the block reader works currently this can only happen if the 'END'
        # is at the very end of the file.
        s = invalid_header("END", False)
        with pytest.warns(
            AstropyUserWarning, match="Missing padding to end of the FITS block"
        ) as w:
            # Don't raise an exception on missing padding, but still produce a
            # warning that the END card is incomplete
            h = fits.Header.fromfile(s, padding=False)
            assert h == horig
        assert len(w) == 1

    def test_invalid_characters(self):
        """
        Test header with invalid characters
        """

        # Generate invalid file with non-ASCII character
        h = fits.Header()
        h["FOO"] = "BAR"
        h["COMMENT"] = "hello"
        hdul = fits.PrimaryHDU(header=h, data=np.arange(5))
        hdul.writeto(self.temp("test.fits"))

        with open(self.temp("test.fits"), "rb") as f:
            out = f.read()
        out = out.replace(b"hello", "héllo".encode("latin1"))
        out = out.replace(b"BAR", "BÀR".encode("latin1"))
        with open(self.temp("test2.fits"), "wb") as f2:
            f2.write(out)

        with pytest.warns(
            AstropyUserWarning,
            match="non-ASCII characters are present in the FITS file",
        ) as w:
            h = fits.getheader(self.temp("test2.fits"))
            assert h["FOO"] == "B?R"
            assert h["COMMENT"] == "h?llo"
        assert len(w) == 1

    def test_unnecessary_move(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/125

        Ensures that a header is not modified when setting the position of a
        keyword that's already in its correct position.
        """

        header = fits.Header([("A", "B"), ("B", "C"), ("C", "D")])

        header.set("B", before=2)
        assert list(header) == ["A", "B", "C"]
        assert not header._modified

        header.set("B", after=0)
        assert list(header) == ["A", "B", "C"]
        assert not header._modified

        header.set("B", before="C")
        assert list(header) == ["A", "B", "C"]
        assert not header._modified

        header.set("B", after="A")
        assert list(header) == ["A", "B", "C"]
        assert not header._modified

        header.set("B", before=2)
        assert list(header) == ["A", "B", "C"]
        assert not header._modified

        # 123 is well past the end, and C is already at the end, so it's in the
        # right place already
        header.set("C", before=123)
        assert list(header) == ["A", "B", "C"]
        assert not header._modified

        header.set("C", after=123)
        assert list(header) == ["A", "B", "C"]
        assert not header._modified

    def test_invalid_float_cards(self):
        """Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/137"""

        # Create a header containing two of the problematic cards in the test
        # case where this came up:
        hstr = "FOCALLEN= +1.550000000000e+002\nAPERTURE= +0.000000000000e+000"
        h = fits.Header.fromstring(hstr, sep="\n")

        # First the case that *does* work prior to fixing this issue
        assert h["FOCALLEN"] == 155.0
        assert h["APERTURE"] == 0.0

        # Now if this were reserialized, would new values for these cards be
        # written with repaired exponent signs?
        with pytest.warns(fits.verify.VerifyWarning) as w:
            assert str(h.cards["FOCALLEN"]) == _pad("FOCALLEN= +1.550000000000E+002")
        assert "Verification reported errors" in str(w[0].message)
        assert h.cards["FOCALLEN"]._modified
        with pytest.warns(fits.verify.VerifyWarning) as w:
            assert str(h.cards["APERTURE"]) == _pad("APERTURE= +0.000000000000E+000")
        assert "Verification reported errors" in str(w[0].message)
        assert h.cards["APERTURE"]._modified
        assert h._modified

        # This is the case that was specifically causing problems; generating
        # the card strings *before* parsing the values.  Also, the card strings
        # really should be "fixed" before being returned to the user
        h = fits.Header.fromstring(hstr, sep="\n")
        with pytest.warns(fits.verify.VerifyWarning) as w:
            assert str(h.cards["FOCALLEN"]) == _pad("FOCALLEN= +1.550000000000E+002")
        assert "Verification reported errors" in str(w[0].message)
        assert h.cards["FOCALLEN"]._modified
        with pytest.warns(fits.verify.VerifyWarning) as w:
            assert str(h.cards["APERTURE"]) == _pad("APERTURE= +0.000000000000E+000")
        assert "Verification reported errors" in str(w[0].message)
        assert h.cards["APERTURE"]._modified

        assert h["FOCALLEN"] == 155.0
        assert h["APERTURE"] == 0.0
        assert h._modified

        # For the heck of it, try assigning the identical values and ensure
        # that the newly fixed value strings are left intact
        h["FOCALLEN"] = 155.0
        h["APERTURE"] = 0.0
        assert str(h.cards["FOCALLEN"]) == _pad("FOCALLEN= +1.550000000000E+002")
        assert str(h.cards["APERTURE"]) == _pad("APERTURE= +0.000000000000E+000")

    def test_invalid_float_cards2(self, capsys):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/140
        """

        # The example for this test requires creating a FITS file containing a
        # slightly misformatted float value.  I can't actually even find a way
        # to do that directly through Astropy--it won't let me.
        hdu = fits.PrimaryHDU()
        hdu.header["TEST"] = 5.0022221e-07
        hdu.writeto(self.temp("test.fits"))

        # Here we manually make the file invalid
        with open(self.temp("test.fits"), "rb+") as f:
            f.seek(346)  # Location of the exponent 'E' symbol
            f.write(encode_ascii("e"))

        with fits.open(self.temp("test.fits")) as hdul, pytest.warns(
            AstropyUserWarning
        ) as w:
            hdul.writeto(self.temp("temp.fits"), output_verify="warn")
        assert len(w) == 5
        # The first two warnings are just the headers to the actual warning
        # message (HDU 0, Card 4).  I'm still not sure things like that
        # should be output as separate warning messages, but that's
        # something to think about...
        msg = str(w[3].message)
        assert "(invalid value string: '5.0022221e-07')" in msg

    def test_leading_zeros(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/137, part 2

        Ticket https://aeon.stsci.edu/ssb/trac/pyfits/ticket/137 also showed that in
        float values like 0.001 the leading zero was unnecessarily being
        stripped off when rewriting the header.  Though leading zeros should be
        removed from integer values to prevent misinterpretation as octal by
        python (for now Astropy will still maintain the leading zeros if now
        changes are made to the value, but will drop them if changes are made).
        """

        c = fits.Card.fromstring("APERTURE= +0.000000000000E+000")
        assert str(c) == _pad("APERTURE= +0.000000000000E+000")
        assert c.value == 0.0
        c = fits.Card.fromstring("APERTURE= 0.000000000000E+000")
        assert str(c) == _pad("APERTURE= 0.000000000000E+000")
        assert c.value == 0.0
        c = fits.Card.fromstring("APERTURE= 017")
        assert str(c) == _pad("APERTURE= 017")
        assert c.value == 17

    def test_assign_boolean(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/123

        Tests assigning Python and Numpy boolean values to keyword values.
        """

        fooimg = _pad("FOO     =                    T")
        barimg = _pad("BAR     =                    F")
        h = fits.Header()
        h["FOO"] = True
        h["BAR"] = False
        assert h["FOO"] is True
        assert h["BAR"] is False
        assert str(h.cards["FOO"]) == fooimg
        assert str(h.cards["BAR"]) == barimg

        h = fits.Header()
        h["FOO"] = np.bool_(True)
        h["BAR"] = np.bool_(False)
        assert h["FOO"] is True
        assert h["BAR"] is False
        assert str(h.cards["FOO"]) == fooimg
        assert str(h.cards["BAR"]) == barimg

        h = fits.Header()
        h.append(fits.Card.fromstring(fooimg))
        h.append(fits.Card.fromstring(barimg))
        assert h["FOO"] is True
        assert h["BAR"] is False
        assert str(h.cards["FOO"]) == fooimg
        assert str(h.cards["BAR"]) == barimg

    def test_header_method_keyword_normalization(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/149

        Basically ensures that all public Header methods are case-insensitive
        w.r.t. keywords.

        Provides a reasonably comprehensive test of several methods at once.
        """

        h = fits.Header([("abC", 1), ("Def", 2), ("GeH", 3)])
        assert list(h) == ["ABC", "DEF", "GEH"]
        assert "abc" in h
        assert "dEf" in h

        assert h["geh"] == 3

        # Case insensitivity of wildcards
        assert len(h["g*"]) == 1

        h["aBc"] = 2
        assert h["abc"] == 2
        # ABC already existed so assigning to aBc should not have added any new
        # cards
        assert len(h) == 3

        del h["gEh"]
        assert list(h) == ["ABC", "DEF"]
        assert len(h) == 2
        assert h.get("def") == 2

        h.set("Abc", 3)
        assert h["ABC"] == 3
        h.set("gEh", 3, before="Abc")
        assert list(h) == ["GEH", "ABC", "DEF"]

        assert h.pop("abC") == 3
        assert len(h) == 2

        assert h.setdefault("def", 3) == 2
        assert len(h) == 2
        assert h.setdefault("aBc", 1) == 1
        assert len(h) == 3
        assert list(h) == ["GEH", "DEF", "ABC"]

        h.update({"GeH": 1, "iJk": 4})
        assert len(h) == 4
        assert list(h) == ["GEH", "DEF", "ABC", "IJK"]
        assert h["GEH"] == 1

        assert h.count("ijk") == 1
        assert h.index("ijk") == 3

        h.remove("Def")
        assert len(h) == 3
        assert list(h) == ["GEH", "ABC", "IJK"]

    def test_end_in_comment(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/142

        Tests a case where the comment of a card ends with END, and is followed
        by several blank cards.
        """

        data = np.arange(100).reshape(10, 10)
        hdu = fits.PrimaryHDU(data=data)
        hdu.header["TESTKW"] = ("Test val", "This is the END")
        # Add a couple blanks after the END string
        hdu.header.append()
        hdu.header.append()
        hdu.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits"), memmap=False) as hdul:
            # memmap = False to avoid leaving open a mmap to the file when we
            # access the data--this causes problems on Windows when we try to
            # overwrite the file later
            assert "TESTKW" in hdul[0].header
            assert hdul[0].header == hdu.header
            assert (hdul[0].data == data).all()

        # Add blanks until the header is extended to two block sizes
        while len(hdu.header) < 36:
            hdu.header.append()
        hdu.writeto(self.temp("test.fits"), overwrite=True)

        with fits.open(self.temp("test.fits")) as hdul:
            assert "TESTKW" in hdul[0].header
            assert hdul[0].header == hdu.header
            assert (hdul[0].data == data).all()

        # Test parsing the same header when it's written to a text file
        hdu.header.totextfile(self.temp("test.hdr"))
        header2 = fits.Header.fromtextfile(self.temp("test.hdr"))
        assert hdu.header == header2

    def test_assign_unicode(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/134

        Assigning a unicode literal as a header value should not fail silently.
        If the value can be converted to ASCII then it should just work.
        Otherwise it should fail with an appropriate value error.

        Also tests unicode for keywords and comments.
        """

        erikku = "\u30a8\u30ea\u30c3\u30af"

        def assign(keyword, val):
            h[keyword] = val

        h = fits.Header()
        h["FOO"] = "BAR"
        assert "FOO" in h
        assert h["FOO"] == "BAR"
        assert repr(h) == _pad("FOO     = 'BAR     '")
        pytest.raises(ValueError, assign, erikku, "BAR")

        h["FOO"] = "BAZ"
        assert h["FOO"] == "BAZ"
        assert repr(h) == _pad("FOO     = 'BAZ     '")
        pytest.raises(ValueError, assign, "FOO", erikku)

        h["FOO"] = ("BAR", "BAZ")
        assert h["FOO"] == "BAR"
        assert h.comments["FOO"] == "BAZ"
        assert repr(h) == _pad("FOO     = 'BAR     '           / BAZ")

        pytest.raises(ValueError, assign, "FOO", ("BAR", erikku))
        pytest.raises(ValueError, assign, "FOO", (erikku, "BAZ"))
        pytest.raises(ValueError, assign, "FOO", (erikku, erikku))

    def test_assign_non_ascii(self):
        """
        First regression test for
        https://github.com/spacetelescope/PyFITS/issues/37

        Although test_assign_unicode ensures that `str` objects containing
        non-ASCII characters cannot be assigned to headers.

        It should not be possible to assign bytes to a header at all.
        """

        h = fits.Header()
        with pytest.raises(ValueError, match="Illegal value: b'Hello'."):
            h.set("TEST", b"Hello")

    def test_header_strip_whitespace(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/146, and
        for the solution that is optional stripping of whitespace from the end
        of a header value.

        By default extra whitespace is stripped off, but if
        `fits.conf.strip_header_whitespace` = False it should not be
        stripped.
        """

        h = fits.Header()
        h["FOO"] = "Bar      "
        assert h["FOO"] == "Bar"
        c = fits.Card.fromstring("QUX     = 'Bar        '")
        h.append(c)
        assert h["QUX"] == "Bar"
        assert h.cards["FOO"].image.rstrip() == "FOO     = 'Bar      '"
        assert h.cards["QUX"].image.rstrip() == "QUX     = 'Bar        '"

        with fits.conf.set_temp("strip_header_whitespace", False):
            assert h["FOO"] == "Bar      "
            assert h["QUX"] == "Bar        "
            assert h.cards["FOO"].image.rstrip() == "FOO     = 'Bar      '"
            assert h.cards["QUX"].image.rstrip() == "QUX     = 'Bar        '"

        assert h["FOO"] == "Bar"
        assert h["QUX"] == "Bar"
        assert h.cards["FOO"].image.rstrip() == "FOO     = 'Bar      '"
        assert h.cards["QUX"].image.rstrip() == "QUX     = 'Bar        '"

    def test_keep_duplicate_history_in_orig_header(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/156

        When creating a new HDU from an existing Header read from an existing
        FITS file, if the original header contains duplicate HISTORY values
        those duplicates should be preserved just as in the original header.

        This bug occurred due to naivete in Header.extend.
        """

        history = [
            "CCD parameters table ...",
            "   reference table oref$n951041ko_ccd.fits",
            "     INFLIGHT 12/07/2001 25/02/2002",
            "     all bias frames",
        ] * 3

        hdu = fits.PrimaryHDU()
        # Add the history entries twice
        for item in history:
            hdu.header["HISTORY"] = item

        hdu.writeto(self.temp("test.fits"))

        with fits.open(self.temp("test.fits")) as hdul:
            assert hdul[0].header["HISTORY"] == history

        new_hdu = fits.PrimaryHDU(header=hdu.header)
        assert new_hdu.header["HISTORY"] == hdu.header["HISTORY"]
        new_hdu.writeto(self.temp("test2.fits"))

        with fits.open(self.temp("test2.fits")) as hdul:
            assert hdul[0].header["HISTORY"] == history

    def test_invalid_keyword_cards(self):
        """
        Test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/109

        Allow opening files with headers containing invalid keywords.
        """

        # Create a header containing a few different types of BAD headers.
        c1 = fits.Card.fromstring("CLFIND2D: contour = 0.30")
        c2 = fits.Card.fromstring("Just some random text.")
        c3 = fits.Card.fromstring("A" * 80)

        hdu = fits.PrimaryHDU()
        # This should work with some warnings
        with pytest.warns(AstropyUserWarning) as w:
            hdu.header.append(c1)
            hdu.header.append(c2)
            hdu.header.append(c3)
        assert len(w) == 3

        hdu.writeto(self.temp("test.fits"))

        with pytest.warns(AstropyUserWarning) as w:
            with fits.open(self.temp("test.fits")) as hdul:
                # Merely opening the file should blast some warnings about the
                # invalid keywords
                assert len(w) == 3

                header = hdul[0].header
                assert "CLFIND2D" in header
                assert "Just som" in header
                assert "AAAAAAAA" in header

                assert header["CLFIND2D"] == ": contour = 0.30"
                assert header["Just som"] == "e random text."
                assert header["AAAAAAAA"] == "A" * 72

                # It should not be possible to assign to the invalid keywords
                pytest.raises(ValueError, header.set, "CLFIND2D", "foo")
                pytest.raises(ValueError, header.set, "Just som", "foo")
                pytest.raises(ValueError, header.set, "AAAAAAAA", "foo")

    def test_fix_hierarch_with_invalid_value(self, capsys):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/172

        Ensures that when fixing a hierarch card it remains a hierarch card.
        """

        c = fits.Card.fromstring("HIERARCH ESO DET CHIP PXSPACE = 5e6")

        with pytest.warns(fits.verify.VerifyWarning) as w:
            c.verify("fix")
        assert "Verification reported errors" in str(w[0].message)
        assert str(c) == _pad("HIERARCH ESO DET CHIP PXSPACE = 5E6")

    def test_assign_inf_nan(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/11

        For the time being it should not be possible to assign the floating
        point values inf or nan to a header value, since this is not defined by
        the FITS standard.
        """

        h = fits.Header()
        pytest.raises(ValueError, h.set, "TEST", float("nan"))
        pytest.raises(ValueError, h.set, "TEST", np.nan)
        pytest.raises(ValueError, h.set, "TEST", np.float32("nan"))
        pytest.raises(ValueError, h.set, "TEST", float("inf"))
        pytest.raises(ValueError, h.set, "TEST", np.inf)

    def test_update_bool(self):
        """
        Regression test for an issue where a value of True in a header
        cannot be updated to a value of 1, and likewise for False/0.
        """

        h = fits.Header([("TEST", True)])
        h["TEST"] = 1
        assert h["TEST"] is not True
        assert isinstance(h["TEST"], int)
        assert h["TEST"] == 1

        h["TEST"] = np.bool_(True)
        assert h["TEST"] is True

        h["TEST"] = False
        assert h["TEST"] is False
        h["TEST"] = np.bool_(False)
        assert h["TEST"] is False

        h["TEST"] = 0
        assert h["TEST"] is not False
        assert isinstance(h["TEST"], int)
        assert h["TEST"] == 0

        h["TEST"] = np.bool_(False)
        assert h["TEST"] is False

    def test_update_numeric(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/49

        Ensure that numeric values can be upcast/downcast between int, float,
        and complex by assigning values that compare equal to the existing
        value but are a different type.
        """

        h = fits.Header()
        h["TEST"] = 1

        # int -> float
        h["TEST"] = 1.0
        assert isinstance(h["TEST"], float)
        assert str(h).startswith("TEST    =                  1.0")

        # float -> int
        h["TEST"] = 1
        assert isinstance(h["TEST"], int)
        assert str(h).startswith("TEST    =                    1")

        # int -> complex
        h["TEST"] = 1.0 + 0.0j
        assert isinstance(h["TEST"], complex)
        assert str(h).startswith("TEST    =           (1.0, 0.0)")

        # complex -> float
        h["TEST"] = 1.0
        assert isinstance(h["TEST"], float)
        assert str(h).startswith("TEST    =                  1.0")

        # float -> complex
        h["TEST"] = 1.0 + 0.0j
        assert isinstance(h["TEST"], complex)
        assert str(h).startswith("TEST    =           (1.0, 0.0)")

        # complex -> int
        h["TEST"] = 1
        assert isinstance(h["TEST"], int)
        assert str(h).startswith("TEST    =                    1")

        # Now the same tests but with zeros
        h["TEST"] = 0

        # int -> float
        h["TEST"] = 0.0
        assert isinstance(h["TEST"], float)
        assert str(h).startswith("TEST    =                  0.0")

        # float -> int
        h["TEST"] = 0
        assert isinstance(h["TEST"], int)
        assert str(h).startswith("TEST    =                    0")

        # int -> complex
        h["TEST"] = 0.0 + 0.0j
        assert isinstance(h["TEST"], complex)
        assert str(h).startswith("TEST    =           (0.0, 0.0)")

        # complex -> float
        h["TEST"] = 0.0
        assert isinstance(h["TEST"], float)
        assert str(h).startswith("TEST    =                  0.0")

        # float -> complex
        h["TEST"] = 0.0 + 0.0j
        assert isinstance(h["TEST"], complex)
        assert str(h).startswith("TEST    =           (0.0, 0.0)")

        # complex -> int
        h["TEST"] = 0
        assert isinstance(h["TEST"], int)
        assert str(h).startswith("TEST    =                    0")

    def test_newlines_in_commentary(self):
        """
        Regression test for https://github.com/spacetelescope/PyFITS/issues/51

        Test data extracted from a header in an actual FITS file found in the
        wild.  Names have been changed to protect the innocent.
        """

        # First ensure that we can't assign new keyword values with newlines in
        # them
        h = fits.Header()
        pytest.raises(ValueError, h.set, "HISTORY", "\n")
        pytest.raises(ValueError, h.set, "HISTORY", "\nabc")
        pytest.raises(ValueError, h.set, "HISTORY", "abc\n")
        pytest.raises(ValueError, h.set, "HISTORY", "abc\ndef")

        test_cards = [
            "HISTORY File modified by user 'wilma' with fv  on 2013-04-22T21:42:18           "
            "HISTORY File modified by user ' fred' with fv  on 2013-04-23T11:16:29           "
            "HISTORY File modified by user ' fred' with fv  on 2013-11-04T16:59:14           "
            "HISTORY File modified by user 'wilma' with fv  on 2013-04-22T21:42:18\nFile modif"
            "HISTORY ied by user 'wilma' with fv  on 2013-04-23T11:16:29\nFile modified by use"
            "HISTORY r ' fred' with fv  on 2013-11-04T16:59:14                               "
            "HISTORY File modified by user 'wilma' with fv  on 2013-04-22T21:42:18\nFile modif"
            "HISTORY ied by user 'wilma' with fv  on 2013-04-23T11:16:29\nFile modified by use"
            "HISTORY r ' fred' with fv  on 2013-11-04T16:59:14\nFile modified by user 'wilma' "
            "HISTORY with fv  on 2013-04-22T21:42:18\nFile modif\nied by user 'wilma' with fv  "
            "HISTORY on 2013-04-23T11:16:29\nFile modified by use\nr ' fred' with fv  on 2013-1"
            "HISTORY 1-04T16:59:14                                                           "
        ]

        for card_image in test_cards:
            c = fits.Card.fromstring(card_image)

            if "\n" in card_image:
                pytest.raises(fits.VerifyError, c.verify, "exception")
            else:
                c.verify("exception")

    def test_long_commentary_card_appended_to_header(self):
        """
        If a HISTORY or COMMENT card with a too-long value is appended to a
        header with Header.append (as opposed to assigning to hdr['HISTORY']
        it fails verification.

        Regression test for https://github.com/astropy/astropy/issues/11486
        """

        header = fits.Header()
        value = "abc" * 90
        # this is what Table does when saving its history metadata key to a
        # FITS file
        header.append(("history", value))
        assert len(header.cards) == 1

        # Test Card._split() directly since this was the main problem area
        key, val = header.cards[0]._split()
        assert key == "HISTORY" and val == value

        # Try writing adding this header to an HDU and writing it to a file
        hdu = fits.PrimaryHDU(header=header)
        hdu.writeto(self.temp("test.fits"), overwrite=True)

    def test_header_fromstring_bytes(self):
        """
        Test reading a Header from a `bytes` string.

        See https://github.com/astropy/astropy/issues/8706
        """

        with open(self.data("test0.fits"), "rb") as fobj:
            pri_hdr_from_bytes = fits.Header.fromstring(fobj.read())

        pri_hdr = fits.getheader(self.data("test0.fits"))
        assert pri_hdr["NAXIS"] == pri_hdr_from_bytes["NAXIS"]
        assert pri_hdr == pri_hdr_from_bytes
        assert pri_hdr.tostring() == pri_hdr_from_bytes.tostring()

    def test_set_keyword_with_space(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/10479
        """
        hdr = fits.Header()
        hdr["KEY2 "] = 2
        hdr["KEY2  "] = 4
        assert len(hdr) == 1
        assert hdr["KEY2"] == 4
        assert hdr["KEY2 "] == 4

    def test_strip(self):
        hdr = fits.getheader(self.data("tb.fits"), ext=1)
        hdr["FOO"] = "bar"
        hdr.strip()
        assert set(hdr) == {"HISTORY", "FOO"}

        hdr = fits.getheader(self.data("tb.fits"), ext=1)
        hdr["FOO"] = "bar"
        hdr = hdr.copy(strip=True)
        assert set(hdr) == {"HISTORY", "FOO"}

    def test_update_invalid_card(self):
        """
        Regression test for https://github.com/astropy/astropy/issues/5408

        Tests updating the value of a card that is malformatted (with an
        invalid value literal).

        This tests two ways of reproducing the problem, one working with a
        Card object directly, and one when reading/writing a header containing
        such an invalid card.
        """

        card = fits.Card.fromstring("KW      = INF  / Comment")
        card.value = "FIXED"
        assert tuple(card) == ("KW", "FIXED", "Comment")
        card.verify("fix")
        assert tuple(card) == ("KW", "FIXED", "Comment")

        card = fits.Card.fromstring("KW      = INF")
        hdu = fits.PrimaryHDU()
        # This is a loophole to write a header containing a malformatted card
        card._verified = True
        hdu.header.append(card)
        hdu.header.tofile(self.temp("bogus.fits"))

        with fits.open(self.temp("bogus.fits")) as hdul:
            hdul[0].header["KW"] = -1
            hdul.writeto(self.temp("bogus_fixed.fits"))

        with fits.open(self.temp("bogus_fixed.fits")) as hdul:
            assert hdul[0].header["KW"] == -1

    def test_index_numpy_int(self):
        header = fits.Header([("A", "FOO"), ("B", 2), ("C", "BAR")])
        idx = np.int8(2)
        assert header[idx] == "BAR"

        header[idx] = "BAZ"
        assert header[idx] == "BAZ"

        header.insert(idx, ("D", 42))
        assert header[idx] == 42

        header.add_comment("HELLO")
        header.add_comment("WORLD")
        assert header["COMMENT"][np.int64(1)] == "WORLD"

        header.append(("C", "BAZBAZ"))
        assert header[("C", np.int16(0))] == "BAZ"
        assert header[("C", np.uint32(1))] == "BAZBAZ"

    def test_header_data_size(self):
        """
        Tests data size calculation (w/o padding) given a Header.
        """
        hdu = fits.PrimaryHDU()
        header = hdu.header
        assert header.data_size == 0

        header["BITPIX"] = 32
        header["NAXIS"] = 2
        header["NAXIS1"] = 100
        header["NAXIS2"] = 100
        assert header.data_size == 40000
        assert header.data_size_padded == 40320


class TestRecordValuedKeywordCards(FitsTestCase):
    """
    Tests for handling of record-valued keyword cards as used by the
    `FITS WCS distortion paper
    <https://www.atnf.csiro.au/people/mcalabre/WCS/dcs_20040422.pdf>`__.

    These tests are derived primarily from the release notes for PyFITS 1.4 (in
    which this feature was first introduced.
    Note that extra leading spaces in the `value` fields should be parsed on input,
    but will be stripped in the cards.
    """

    def setup_method(self):
        super().setup_method()
        self._test_header = fits.Header()
        self._test_header.set("DP1", "NAXIS: 2")
        self._test_header.set("DP1", "AXIS.1: 1")
        self._test_header.set("DP1", "AXIS.2: 2")
        self._test_header.set("DP1", "NAUX:   2")
        self._test_header.set("DP1", "AUX.1.COEFF.0: 0")
        self._test_header.set("DP1", "AUX.1.POWER.0: 1")
        self._test_header.set("DP1", "AUX.1.COEFF.1: 0.00048828125")
        self._test_header.set("DP1", "AUX.1.POWER.1:  1")

    def test_initialize_rvkc(self):
        """
        Test different methods for initializing a card that should be
        recognized as a RVKC
        """

        c = fits.Card.fromstring("DP1     = 'NAXIS: 2' / A comment")
        assert c.keyword == "DP1.NAXIS"
        assert c.value == 2.0
        assert c.field_specifier == "NAXIS"
        assert c.comment == "A comment"

        c = fits.Card.fromstring("DP1     = 'NAXIS:  2.1'")
        assert c.keyword == "DP1.NAXIS"
        assert c.value == 2.1
        assert c.field_specifier == "NAXIS"

        c = fits.Card.fromstring("DP1     = 'NAXIS: a'")
        assert c.keyword == "DP1"
        assert c.value == "NAXIS: a"
        assert c.field_specifier is None

        c = fits.Card("DP1", "NAXIS: 2")
        assert c.keyword == "DP1.NAXIS"
        assert c.value == 2.0
        assert c.field_specifier == "NAXIS"

        c = fits.Card("DP1", "NAXIS:  2.0")
        assert c.keyword == "DP1.NAXIS"
        assert c.value == 2.0
        assert c.field_specifier == "NAXIS"

        c = fits.Card("DP1", "NAXIS: a")
        assert c.keyword == "DP1"
        assert c.value == "NAXIS: a"
        assert c.field_specifier is None

        c = fits.Card("DP1.NAXIS", 2)
        assert c.keyword == "DP1.NAXIS"
        assert c.value == 2.0
        assert c.field_specifier == "NAXIS"

        c = fits.Card("DP1.NAXIS", 2.0)
        assert c.keyword == "DP1.NAXIS"
        assert c.value == 2.0
        assert c.field_specifier == "NAXIS"

        with pytest.warns(fits.verify.VerifyWarning):
            c = fits.Card("DP1.NAXIS", "a")
        assert c.keyword == "DP1.NAXIS"
        assert c.value == "a"
        assert c.field_specifier is None

    def test_parse_field_specifier(self):
        """
        Tests that the field_specifier can accessed from a card read from a
        string before any other attributes are accessed.
        """

        c = fits.Card.fromstring("DP1     = 'NAXIS: 2' / A comment")
        assert c.field_specifier == "NAXIS"
        assert c.keyword == "DP1.NAXIS"
        assert c.value == 2.0
        assert c.comment == "A comment"

    def test_update_field_specifier(self):
        """
        Test setting the field_specifier attribute and updating the card image
        to reflect the new value.
        """

        c = fits.Card.fromstring("DP1     = 'NAXIS: 2' / A comment")
        assert c.field_specifier == "NAXIS"
        c.field_specifier = "NAXIS1"
        assert c.field_specifier == "NAXIS1"
        assert c.keyword == "DP1.NAXIS1"
        assert c.value == 2.0
        assert c.comment == "A comment"
        assert str(c).rstrip() == "DP1     = 'NAXIS1: 2' / A comment"

    def test_field_specifier_case_senstivity(self):
        """
        The keyword portion of an RVKC should still be case-insensitive, but
        the field-specifier portion should be case-sensitive.
        """

        header = fits.Header()
        header.set("abc.def", 1)
        header.set("abc.DEF", 2)
        assert header["abc.def"] == 1
        assert header["ABC.def"] == 1
        assert header["aBc.def"] == 1
        assert header["ABC.DEF"] == 2
        assert "ABC.dEf" not in header

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
        Returning a RVKC just via the keyword name should return the full value
        string of the first card with that keyword.

        This test was changed to reflect the requirement in ticket
        https://aeon.stsci.edu/ssb/trac/pyfits/ticket/184--previously it required
        _test_header['DP1'] to return the parsed float value.
        """

        assert self._test_header["DP1"] == "NAXIS: 2"

    def test_get_rvkc_by_keyword_and_field_specifier(self):
        """
        Returning a RVKC via the full keyword/field-specifier combination
        should return the floating point value associated with the RVKC.
        """

        assert self._test_header["DP1.NAXIS"] == 2.0
        assert isinstance(self._test_header["DP1.NAXIS"], float)
        assert self._test_header["DP1.AUX.1.COEFF.1"] == 0.00048828125

    def test_access_nonexistent_rvkc(self):
        """
        Accessing a nonexistent RVKC should raise an IndexError for
        index-based lookup, or a KeyError for keyword lookup (like a normal
        card).
        """

        pytest.raises(IndexError, lambda x: self._test_header[x], 8)
        # Test exception with message
        with pytest.raises(KeyError, match=r"Keyword 'DP1\.AXIS\.3' not found."):
            self._test_header["DP1.AXIS.3"]

    def test_update_rvkc(self):
        """A RVKC can be updated either via index or keyword access."""

        self._test_header[0] = 3
        assert self._test_header["DP1.NAXIS"] == 3.0
        assert isinstance(self._test_header["DP1.NAXIS"], float)

        self._test_header["DP1.AXIS.1"] = 1.1
        assert self._test_header["DP1.AXIS.1"] == 1.1

    def test_update_rvkc_2(self):
        """Regression test for an issue that appeared after SVN r2412."""

        h = fits.Header()
        h["D2IM1.EXTVER"] = 1
        assert h["D2IM1.EXTVER"] == 1.0
        h["D2IM1.EXTVER"] = 2
        assert h["D2IM1.EXTVER"] == 2.0

    def test_raw_keyword_value(self):
        c = fits.Card.fromstring("DP1     = 'NAXIS: 2' / A comment")
        assert c.rawkeyword == "DP1"
        assert c.rawvalue == "NAXIS: 2"

        c = fits.Card("DP1.NAXIS", 2)
        assert c.rawkeyword == "DP1"
        assert c.rawvalue == "NAXIS: 2.0"

        c = fits.Card("DP1.NAXIS", 2.0)
        assert c.rawkeyword == "DP1"
        assert c.rawvalue == "NAXIS: 2.0"

    def test_rvkc_insert_after(self):
        """
        It should be possible to insert a new RVKC after an existing one
        specified by the full keyword/field-specifier combination."""

        self._test_header.set("DP1", "AXIS.3: 1", "a comment", after="DP1.AXIS.2")
        assert self._test_header[3] == 1
        assert self._test_header["DP1.AXIS.3"] == 1

    def test_rvkc_delete(self):
        """
        Deleting a RVKC should work as with a normal card by using the full
        keyword/field-spcifier combination.
        """

        del self._test_header["DP1.AXIS.1"]
        assert len(self._test_header) == 7
        assert list(self._test_header)[0] == "DP1.NAXIS"
        assert self._test_header[0] == 2
        assert list(self._test_header)[1] == "DP1.AXIS.2"

        # Perform a subsequent delete to make sure all the index mappings were
        # updated
        del self._test_header["DP1.AXIS.2"]
        assert len(self._test_header) == 6
        assert list(self._test_header)[0] == "DP1.NAXIS"
        assert self._test_header[0] == 2
        assert list(self._test_header)[1] == "DP1.NAUX"
        assert self._test_header[1] == 2

    def test_pattern_matching_keys(self):
        """Test the keyword filter strings with RVKCs."""

        cl = self._test_header["DP1.AXIS.*"]
        assert isinstance(cl, fits.Header)
        assert [str(c).strip() for c in cl.cards] == [
            "DP1     = 'AXIS.1: 1'",
            "DP1     = 'AXIS.2: 2'",
        ]

        cl = self._test_header["DP1.N*"]
        assert [str(c).strip() for c in cl.cards] == [
            "DP1     = 'NAXIS: 2'",
            "DP1     = 'NAUX: 2'",
        ]

        cl = self._test_header["DP1.AUX..."]
        assert [str(c).strip() for c in cl.cards] == [
            "DP1     = 'AUX.1.COEFF.0: 0'",
            "DP1     = 'AUX.1.POWER.0: 1'",
            "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
            "DP1     = 'AUX.1.POWER.1: 1'",
        ]

        cl = self._test_header["DP?.NAXIS"]
        assert [str(c).strip() for c in cl.cards] == ["DP1     = 'NAXIS: 2'"]

        cl = self._test_header["DP1.A*S.*"]
        assert [str(c).strip() for c in cl.cards] == [
            "DP1     = 'AXIS.1: 1'",
            "DP1     = 'AXIS.2: 2'",
        ]

    def test_pattern_matching_key_deletion(self):
        """Deletion by filter strings should work."""

        del self._test_header["DP1.A*..."]
        assert len(self._test_header) == 2
        assert list(self._test_header)[0] == "DP1.NAXIS"
        assert self._test_header[0] == 2
        assert list(self._test_header)[1] == "DP1.NAUX"
        assert self._test_header[1] == 2

    def test_successive_pattern_matching(self):
        """
        A card list returned via a filter string should be further filterable.
        """

        cl = self._test_header["DP1.A*..."]
        assert [str(c).strip() for c in cl.cards] == [
            "DP1     = 'AXIS.1: 1'",
            "DP1     = 'AXIS.2: 2'",
            "DP1     = 'AUX.1.COEFF.0: 0'",
            "DP1     = 'AUX.1.POWER.0: 1'",
            "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
            "DP1     = 'AUX.1.POWER.1: 1'",
        ]

        cl2 = cl["*.*AUX..."]
        assert [str(c).strip() for c in cl2.cards] == [
            "DP1     = 'AUX.1.COEFF.0: 0'",
            "DP1     = 'AUX.1.POWER.0: 1'",
            "DP1     = 'AUX.1.COEFF.1: 0.00048828125'",
            "DP1     = 'AUX.1.POWER.1: 1'",
        ]

    def test_rvkc_in_cardlist_keys(self):
        """
        The CardList.keys() method should return full keyword/field-spec values
        for RVKCs.
        """

        cl = self._test_header["DP1.AXIS.*"]
        assert list(cl) == ["DP1.AXIS.1", "DP1.AXIS.2"]

    def test_rvkc_in_cardlist_values(self):
        """
        The CardList.values() method should return the values of all RVKCs as
        floating point values.
        """

        cl = self._test_header["DP1.AXIS.*"]
        assert list(cl.values()) == [1.0, 2.0]

    def test_rvkc_value_attribute(self):
        """
        Individual card values should be accessible by the .value attribute
        (which should return a float).
        """

        cl = self._test_header["DP1.AXIS.*"]
        assert cl.cards[0].value == 1.0
        assert isinstance(cl.cards[0].value, float)

    def test_overly_permissive_parsing(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/183

        Ensures that cards with standard commentary keywords are never treated
        as RVKCs.  Also ensures that cards not strictly matching the RVKC
        pattern are not treated as such.
        """

        h = fits.Header()
        h["HISTORY"] = "AXIS.1: 2"
        h["HISTORY"] = "AXIS.2: 2"
        assert "HISTORY.AXIS" not in h
        assert "HISTORY.AXIS.1" not in h
        assert "HISTORY.AXIS.2" not in h
        assert h["HISTORY"] == ["AXIS.1: 2", "AXIS.2: 2"]

        # This is an example straight out of the ticket where everything after
        # the '2012' in the date value was being ignored, allowing the value to
        # successfully be parsed as a "float"
        h = fits.Header()
        h["HISTORY"] = "Date: 2012-09-19T13:58:53.756061"
        assert "HISTORY.Date" not in h
        assert str(h.cards[0]) == _pad("HISTORY Date: 2012-09-19T13:58:53.756061")

        c = fits.Card.fromstring("        'Date: 2012-09-19T13:58:53.756061'")
        assert c.keyword == ""
        assert c.value == "'Date: 2012-09-19T13:58:53.756061'"
        assert c.field_specifier is None

        h = fits.Header()
        h["FOO"] = "Date: 2012-09-19T13:58:53.756061"
        assert "FOO.Date" not in h
        assert str(h.cards[0]) == _pad("FOO     = 'Date: 2012-09-19T13:58:53.756061'")

    def test_overly_aggressive_rvkc_lookup(self):
        """
        Regression test for https://aeon.stsci.edu/ssb/trac/pyfits/ticket/184

        Ensures that looking up a RVKC by keyword only (without the
        field-specifier) in a header returns the full string value of that card
        without parsing it as a RVKC.  Also ensures that a full field-specifier
        is required to match a RVKC--a partial field-specifier that doesn't
        explicitly match any record-valued keyword should result in a KeyError.
        """

        c1 = fits.Card.fromstring("FOO     = 'AXIS.1: 2'")
        c2 = fits.Card.fromstring("FOO     = 'AXIS.2: 4'")
        h = fits.Header([c1, c2])
        assert h["FOO"] == "AXIS.1: 2"
        assert h[("FOO", 1)] == "AXIS.2: 4"
        assert h["FOO.AXIS.1"] == 2.0
        assert h["FOO.AXIS.2"] == 4.0
        assert "FOO.AXIS" not in h
        assert "FOO.AXIS." not in h
        assert "FOO." not in h
        pytest.raises(KeyError, lambda: h["FOO.AXIS"])
        pytest.raises(KeyError, lambda: h["FOO.AXIS."])
        pytest.raises(KeyError, lambda: h["FOO."])

    def test_fitsheader_script(self):
        """Tests the basic functionality of the `fitsheader` script."""
        from astropy.io.fits.scripts import fitsheader

        # Can an extension by specified by the EXTNAME keyword?
        hf = fitsheader.HeaderFormatter(self.data("zerowidth.fits"))
        output = hf.parse(extensions=["AIPS FQ"])
        assert "EXTNAME = 'AIPS FQ" in output
        assert "BITPIX" in output

        # Can we limit the display to one specific keyword?
        output = hf.parse(extensions=["AIPS FQ"], keywords=["EXTNAME"])
        assert "EXTNAME = 'AIPS FQ" in output
        assert "BITPIX  =" not in output
        assert len(output.split("\n")) == 3

        # Can we limit the display to two specific keywords?
        output = hf.parse(extensions=[1], keywords=["EXTNAME", "BITPIX"])
        assert "EXTNAME =" in output
        assert "BITPIX  =" in output
        assert len(output.split("\n")) == 4

        # Can we use wildcards for keywords?
        output = hf.parse(extensions=[1], keywords=["NAXIS*"])
        assert "NAXIS   =" in output
        assert "NAXIS1  =" in output
        assert "NAXIS2  =" in output
        hf.close()

        # Can an extension by specified by the EXTNAME+EXTVER keywords?
        hf = fitsheader.HeaderFormatter(self.data("test0.fits"))
        assert "EXTNAME = 'SCI" in hf.parse(extensions=["SCI,2"])
        hf.close()

        # Can we print the original header before decompression?
        hf = fitsheader.HeaderFormatter(self.data("comp.fits"))
        assert "XTENSION= 'IMAGE" in hf.parse(extensions=[1], compressed=False)
        assert "XTENSION= 'BINTABLE" in hf.parse(extensions=[1], compressed=True)
        hf.close()

    def test_fitsheader_compressed_from_primary_image_ext(self):
        """Regression test for issue https://github.com/astropy/astropy/issues/7312"""
        data = np.arange(2 * 2, dtype=np.int8).reshape((2, 2))
        phdu = fits.PrimaryHDU(data=data)
        chdu = fits.CompImageHDU(data=phdu.data, header=phdu.header)
        chdu.writeto(self.temp("tmp2.fits"), overwrite=True)
        with fits.open(self.temp("tmp2.fits")) as hdul:
            assert "XTENSION" not in hdul[1].header
            assert "PCOUNT" not in hdul[1].header
            assert "GCOUNT" not in hdul[1].header

    def test_fitsheader_table_feature(self):
        """Tests the `--table` feature of the `fitsheader` script."""
        from astropy.io import fits
        from astropy.io.fits.scripts import fitsheader

        test_filename = self.data("zerowidth.fits")

        formatter = fitsheader.TableHeaderFormatter(test_filename)

        with fits.open(test_filename) as fitsobj:
            # Does the table contain the expected number of rows?
            mytable = formatter.parse([0])
            assert len(mytable) == len(fitsobj[0].header)
            # Repeat the above test when multiple HDUs are requested
            mytable = formatter.parse(extensions=["AIPS FQ", 2, "4"])
            assert len(mytable) == (
                len(fitsobj["AIPS FQ"].header)
                + len(fitsobj[2].header)
                + len(fitsobj[4].header)
            )

        # Can we recover the filename and extension name from the table?
        mytable = formatter.parse(extensions=["AIPS FQ"])
        assert np.all(mytable["filename"] == test_filename)
        assert np.all(mytable["hdu"] == "AIPS FQ")
        assert mytable["value"][mytable["keyword"] == "EXTNAME"] == "AIPS FQ"

        # Can we specify a single extension/keyword?
        mytable = formatter.parse(extensions=["AIPS FQ"], keywords=["EXTNAME"])
        assert len(mytable) == 1
        assert mytable["hdu"][0] == "AIPS FQ"
        assert mytable["keyword"][0] == "EXTNAME"
        assert mytable["value"][0] == "AIPS FQ"

        # Is an incorrect extension dealt with gracefully?
        mytable = formatter.parse(extensions=["DOES_NOT_EXIST"])
        assert mytable is None
        # Is an incorrect keyword dealt with gracefully?
        mytable = formatter.parse(extensions=["AIPS FQ"], keywords=["DOES_NOT_EXIST"])
        assert mytable is None
        formatter.close()

    @pytest.mark.parametrize("mode", ["wb", "wb+", "ab", "ab+"])
    def test_hdu_writeto_mode(self, mode):
        with open(self.temp("mode.fits"), mode=mode) as ff:
            hdu = fits.ImageHDU(data=np.ones(5))
            hdu.writeto(ff)


def test_subclass():
    """Check that subclasses don't get ignored on slicing and copying."""

    class MyHeader(fits.Header):
        def append(self, card, *args, **kwargs):
            if isinstance(card, tuple) and len(card) == 2:
                # Just for our checks we add a comment if there is none.
                card += ("no comment",)

            return super().append(card, *args, **kwargs)

    my_header = MyHeader(
        (
            ("a", 1.0, "first"),
            ("b", 2.0, "second"),
            (
                "c",
                3.0,
            ),
        )
    )

    assert my_header.comments["a"] == "first"
    assert my_header.comments["b"] == "second"
    assert my_header.comments["c"] == "no comment"

    slice_ = my_header[1:]
    assert type(slice_) is MyHeader
    assert slice_.comments["b"] == "second"
    assert slice_.comments["c"] == "no comment"
    selection = my_header["c*"]
    assert type(selection) is MyHeader
    assert selection.comments["c"] == "no comment"
    copy_ = my_header.copy()
    assert type(copy_) is MyHeader
    assert copy_.comments["b"] == "second"
    assert copy_.comments["c"] == "no comment"
    my_header.extend((("d", 4.0),))
    assert my_header.comments["d"] == "no comment"
