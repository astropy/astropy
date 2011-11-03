from __future__ import division # confidence high
from __future__ import with_statement

import os
import warnings

import numpy as np

import pyfits
from pyfits.tests import PyfitsTestCase
from pyfits.tests.util import CaptureStdout, catch_warnings

from nose.tools import assert_equal, assert_raises, assert_true, assert_false


class TestImageFunctions(PyfitsTestCase):
    def test_card_constructor_default_args(self):
        """Test the constructor with default argument values."""

        c = pyfits.Card()
        assert_equal('', c.key)

    def test_constructor_name_arg(self):
        """Like the test of the same name in test_table.py"""

        hdu = pyfits.ImageHDU()
        assert_equal(hdu.name, '')
        assert_true('EXTNAME' not in hdu.header)
        hdu.name = 'FOO'
        assert_equal(hdu.name, 'FOO')
        assert_equal(hdu.header['EXTNAME'], 'FOO')

        # Passing name to constructor
        hdu = pyfits.ImageHDU(name='FOO')
        assert_equal(hdu.name, 'FOO')
        assert_equal(hdu.header['EXTNAME'], 'FOO')

        # And overriding a header with a different extname
        hdr = pyfits.Header()
        hdr.update('EXTNAME', 'EVENTS')
        hdu = pyfits.ImageHDU(header=hdr, name='FOO')
        assert_equal(hdu.name, 'FOO')
        assert_equal(hdu.header['EXTNAME'], 'FOO')

    def test_fromstring_set_attribute_ascardimage(self):
        """Test fromstring() which will return a new card."""

        c = pyfits.Card('abc', 99).fromstring('xyz     = 100')
        assert_equal(100, c.value)

        # test set attribute and  ascardimage() using the most updated attributes
        c.value = 200
        assert_equal(c.ascardimage(),
                     "XYZ     =                  200                                                  ")

    def test_string(self):
        """Test string value"""

        c = pyfits.Card('abc', '<8 ch')
        assert_equal(str(c),
                     "ABC     = '<8 ch   '                                                            ")
        c = pyfits.Card('nullstr', '')
        assert_equal(str(c),
                     "NULLSTR = ''                                                                    ")

    def test_boolean_value_card(self):
        """Boolean value card"""

        c = pyfits.Card("abc", True)
        assert_equal(str(c),
                     "ABC     =                    T                                                  ")

        c = pyfits.Card.fromstring('abc     = F')
        assert_equal(c.value, False)

    def test_long_integer_number(self):
        """long integer number"""

        c = pyfits.Card('long_int', -467374636747637647347374734737437)
        assert_equal(str(c),
                     "LONG_INT= -467374636747637647347374734737437                                    ")

    def test_floating_point_number(self):
        """ floating point number"""

        c = pyfits.Card('floatnum', -467374636747637647347374734737437.)

        if (str(c) != "FLOATNUM= -4.6737463674763E+32                                                  " and
            str(c) != "FLOATNUM= -4.6737463674763E+032                                                 "):
            assert_equal(str(c),
                         "FLOATNUM= -4.6737463674763E+32                                                  ")

    def test_complex_value(self):
        """complex value"""

        c = pyfits.Card('abc',
                        1.2345377437887837487e88+6324767364763746367e-33j)

        if (str(c) != "ABC     = (1.23453774378878E+88, 6.32476736476374E-15)                          " and
            str(c) != "ABC     = (1.2345377437887E+088, 6.3247673647637E-015)                          "):
            assert_equal(str(c),
                         "ABC     = (1.23453774378878E+88, 6.32476736476374E-15)                          ")

    def test_card_image_constructed_too_long(self):
        with CaptureStdout():
            # card image constructed from key/value/comment is too long
            # (non-string value)
            c = pyfits.Card('abc', 9, 'abcde'*20)
            assert_equal(str(c),
                         "ABC     =                    9 / abcdeabcdeabcdeabcdeabcdeabcdeabcdeabcdeabcdeab")
            c = pyfits.Card('abc', 'a'*68, 'abcdefg')
            assert_equal(str(c),
                         "ABC     = 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'")

    def test_constructor_filter_illegal_data_structures(self):
        # the constrctor will filter out illegal data structures...
        assert_raises(ValueError, pyfits.Card, ('abc',), {'value': (2, 3)})
        assert_raises(ValueError, pyfits.Card, 'key', [], 'comment')

    def test_keyword_too_long(self):
        #keywords too long
        assert_raises(ValueError, pyfits.Card, 'abcdefghi', 'long')


    def test_illegal_characters_in_key(self):
        # will not allow illegal characters in key when using constructor
        assert_raises(ValueError, pyfits.Card, 'abc+', 9)


    def test_ascardiage_verifies_the_comment_string_to_be_ascii_text(self):
        # the ascardimage() verifies the comment string to be ASCII text
        c = pyfits.Card.fromstring('abc     = +  2.1   e + 12 / abcde\0')
        assert_raises(Exception, c.ascardimage)

    def test_commentary_cards(self):
        # commentary cards
        c = pyfits.Card("history",
                        "A commentary card's value has no quotes around it.")
        assert_equal(str(c),
                     "HISTORY A commentary card's value has no quotes around it.                      ")
        c = pyfits.Card("comment",
                        "A commentary card has no comment.", "comment")
        assert_equal(str(c),
                     "COMMENT A commentary card has no comment.                                       ")

    def test_commentary_card_created_by_fromstring(self):
        # commentary card created by fromstring()
        c = pyfits.Card.fromstring("COMMENT card has no comments. / text after slash is still part of the value.")
        assert_equal(c.value,
                     'card has no comments. / text after slash is still part of the value.')
        assert_equal(c.comment, '')

    def test_commentary_card_will_not_parse_numerical_value(self):
        # commentary card will not parse the numerical value
        c = pyfits.Card.fromstring("history  (1, 2)")
        assert_equal(str(c.ascardimage()),
                     "HISTORY  (1, 2)                                                                 ")

    def test_equal_sign_after_column8(self):
        # equal sign after column 8 of a commentary card will be part ofthe string value
        c = pyfits.Card.fromstring("history =   (1, 2)")
        assert_equal(str(c.ascardimage()),
                     "HISTORY =   (1, 2)                                                              ")

    def test_specify_undefined_value(self):
        # this is how to specify an undefined value
        c = pyfits.Card("undef", pyfits.card.UNDEFINED)
        assert_equal(str(c),
                     "UNDEF   =                                                                       ")

    def test_complex_number_using_string_input(self):
        # complex number using string input
        c = pyfits.Card.fromstring('abc     = (8, 9)')
        assert_equal(str(c.ascardimage()),
                     "ABC     =               (8, 9)                                                  ")

    def test_fixable_non_standard_fits_card(self):
        # fixable non-standard FITS card will keep the original format
        c = pyfits.Card.fromstring('abc     = +  2.1   e + 12')
        assert_equal(c.value,2100000000000.0)
        assert_equal(str(c.ascardimage()),
                     "ABC     =             +2.1E+12                                                  ")

    def test_fixable_non_fsc(self):
        # fixable non-FSC: if the card is not parsable, it's value will be
        # assumed
        # to be a string and everything after the first slash will be comment
        c = pyfits.Card.fromstring("no_quote=  this card's value has no quotes / let's also try the comment")
        assert_equal(str(c.ascardimage()),
                     "NO_QUOTE= 'this card''s value has no quotes' / let's also try the comment       ")

    def test_undefined_value_using_string_input(self):
        # undefined value using string input
        c = pyfits.Card.fromstring('abc     =    ')
        assert_equal(str(c.ascardimage()),
                     "ABC     =                                                                       ")

    def test_misalocated_equal_sign(self):
        # test mislocated "=" sign
        c = pyfits.Card.fromstring('xyz= 100')
        assert_equal(c.key, 'XYZ')
        assert_equal(c.value, 100)
        assert_equal(str(c.ascardimage()),
                     "XYZ     =                  100                                                  ")

    def test_equal_only_up_to_column_10(self):
        # the test of "=" location is only up to column 10
        c = pyfits.Card.fromstring("histo       =   (1, 2)")
        assert_equal(str(c.ascardimage()),
                     "HISTO   = '=   (1, 2)'                                                          ")
        c = pyfits.Card.fromstring("   history          (1, 2)")
        assert_equal(str(c.ascardimage()),
                     "HISTO   = 'ry          (1, 2)'                                                  ")

    def test_verification(self):
        # verification
        c = pyfits.Card.fromstring('abc= a6')
        with CaptureStdout() as f:
            c.verify()
            assert_true(
                'Card image is not FITS standard (equal sign not at column 8).'
                in f.getvalue())
        assert_equal(str(c),
                     "abc= a6                                                                         ")

    def test_fix(self):
        c = pyfits.Card.fromstring('abc= a6')
        with CaptureStdout() as f:
            c.verify('fix')
            assert_true('Fixed card to be FITS standard.: ABC' in f.getvalue())
        assert_equal(str(c),
                     "ABC     = 'a6      '                                                            ")

    def test_long_string_value(self):
        # test long string value
        c = pyfits.Card('abc', 'long string value '*10, 'long comment '*10)
        assert_equal(str(c),
            "ABC     = 'long string value long string value long string value long string &' "
            "CONTINUE  'value long string value long string value long string value long &'  "
            "CONTINUE  'string value long string value long string value &'                  "
            "CONTINUE  '&' / long comment long comment long comment long comment long        "
            "CONTINUE  '&' / comment long comment long comment long comment long comment     "
            "CONTINUE  '&' / long comment                                                    ")

    def test_long_string_from_file(self):
        c = pyfits.Card('abc', 'long string value '*10, 'long comment '*10)
        hdu = pyfits.PrimaryHDU()
        hdu.header.ascard.append(c)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = pyfits.open(self.temp('test_new.fits'))
        c = hdul[0].header.ascard['abc']
        hdul.close()
        assert_equal(str(c),
            "ABC     = 'long string value long string value long string value long string &' "
            "CONTINUE  'value long string value long string value long string value long &'  "
            "CONTINUE  'string value long string value long string value &'                  "
            "CONTINUE  '&' / long comment long comment long comment long comment long        "
            "CONTINUE  '&' / comment long comment long comment long comment long comment     "
            "CONTINUE  '&' / long comment                                                    ")


    def test_word_in_long_string_too_long(self):
        # if a word in a long string is too long, it will be cut in the middle
        c = pyfits.Card('abc', 'longstringvalue'*10, 'longcomment'*10)
        assert_equal(str(c),
            "ABC     = 'longstringvaluelongstringvaluelongstringvaluelongstringvaluelongstr&'"
            "CONTINUE  'ingvaluelongstringvaluelongstringvaluelongstringvaluelongstringvalu&'"
            "CONTINUE  'elongstringvalue&'                                                   "
            "CONTINUE  '&' / longcommentlongcommentlongcommentlongcommentlongcommentlongcomme"
            "CONTINUE  '&' / ntlongcommentlongcommentlongcommentlongcomment                  ")

    def test_long_string_value_via_fromstring(self):
        # long string value via fromstring() method
        c = pyfits.Card.fromstring(
            pyfits.card._pad("abc     = 'longstring''s testing  &  ' / comments in line 1") +
            pyfits.card._pad("continue  'continue with long string but without the ampersand at the end' /") +
            pyfits.card._pad("continue  'continue must have string value (with quotes)' / comments with ''. "))
        assert_equal(str(c.ascardimage()),
            "ABC     = 'longstring''s testing  continue with long string but without the &'  "
            "CONTINUE  'ampersand at the endcontinue must have string value (with quotes)&'  "
            "CONTINUE  '&' / comments in line 1 comments with ''.                            ")

    def test_hierarch_card(self):
        c = pyfits.Card('hierarch abcdefghi', 10)
        assert_equal(str(c.ascardimage()),
            "HIERARCH abcdefghi = 10                                                         ")
        c = pyfits.Card('HIERARCH ESO INS SLIT2 Y1FRML', 'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)')
        assert_equal(str(c.ascardimage()),
            "HIERARCH ESO INS SLIT2 Y1FRML= 'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)'")


    def test_open(self):
        # The function "open" reads a FITS file into an HDUList object.  There
        # are three modes to open: "readonly" (the default), "append", and
        # "update".

        # Open a file read-only (the default mode), the content of the FITS
        # file are read into memory.
        r = pyfits.open(self.data('test0.fits')) # readonly

        # data parts are latent instantiation, so if we close the HDUList
        # without touching data, data can not be accessed.
        r.close()
        assert_raises(Exception, lambda x: x[1].data[:2,:2], r)

    def test_open_2(self):
        r = pyfits.open(self.data('test0.fits'))

        info = [(0, 'PRIMARY', 'PrimaryHDU', 138, (), 'int16', '')] + \
               [(x, 'SCI', 'ImageHDU', 61, (40, 40), 'int16', '')
                for x in range(1, 5)]

        try:
            assert_equal(r.info(output=False), info)
        finally:
            r.close()

    def test_io_manipulation(self):
        # Get a keyword value.  An extension can be referred by name or by
        # number.  Both extension and keyword names are case insensitive.
        r = pyfits.open(self.data('test0.fits'))
        assert_equal(r['primary'].header['naxis'], 0)
        assert_equal(r[0].header['naxis'], 0)

        # If there are more than one extension with the same EXTNAME value, the
        # EXTVER can be used (as the second argument) to distinguish the
        # extension.
        assert_equal(r['sci',1].header['detector'], 1)

        # append (using "update()") a new card
        r[0].header.update('xxx', 1.234e56)

        if str(r[0].header.ascard[-3:]) != \
           "EXPFLAG = 'NORMAL            ' / Exposure interruption indicator                \n" \
           "FILENAME= 'vtest3.fits'        / File name                                      \n" \
           "XXX     =            1.234E+56                                                  " and \
           str(r[0].header.ascard[-3:]) != \
           "EXPFLAG = 'NORMAL            ' / Exposure interruption indicator                \n" \
           "FILENAME= 'vtest3.fits'        / File name                                      \n" \
           "XXX     =           1.234E+056                                                  ":
            assert_equal(str(r[0].header.ascard[-3:]),
                "EXPFLAG = 'NORMAL            ' / Exposure interruption indicator                \n"
                "FILENAME= 'vtest3.fits'        / File name                                      \n"
                "XXX     =            1.234E+56                                                  ")

        # rename a keyword
        r[0].header.rename_key('filename', 'fname')
        assert_raises(ValueError, r[0].header.rename_key,
                          'fname', 'history')

        assert_raises(ValueError, r[0].header.rename_key,
                          'fname', 'simple')
        r[0].header.rename_key('fname', 'filename')

        # get a subsection of data
        assert_equal(r[2].data[:3,:3].all(),
                         np.array([[349, 349, 348],
                                   [348, 348, 348],
                                   [349, 349, 350]], dtype=np.int16).all())

        # We can create a new FITS file by opening a new file with "append"
        # mode.
        n=pyfits.open(self.temp('test_new.fits'), mode='append')

        # Append the primary header and the 2nd extension to the new file.
        n.append(r[0])
        n.append(r[2])

        # The flush method will write the current HDUList object back to the
        # newly created file on disk.  The HDUList is still open and can be
        # further operated.
        n.flush()
        assert_equal(n[1].data[1,1], 349)

        #modify a data point
        n[1].data[1,1] = 99

        # When the file is closed, the most recent additions of extension(s)
        # since last flush() will be appended, but any HDU already existed at
        # the last flush will not be modified
        n.close()

        # If an existing file is opened with "append" mode, like the readonly
        # mode, the HDU's will be read into the HDUList which can be modified
        # in memory but can not be written back to the original file.  A file
        # opened with append mode can only add new HDU's.
        os.rename(self.temp('test_new.fits'), self.temp('test_append.fits'))
        a = pyfits.open(self.temp('test_append.fits'), mode='append')

        # The above change did not take effect since this was made after the
        # flush().
        assert_equal(a[1].data[1,1], 349)

        a.append(r[1])
        a.close()

        # When changes are made to an HDUList which was opened with "update"
        # mode, they will be written back to the original file when a
        # flush/close is called.
        os.rename(self.temp('test_append.fits'), self.temp('test_update.fits'))

        u = pyfits.open(self.temp('test_update.fits'), mode='update')

        # When the changes do not alter the size structures of the original (or
        # since last flush) HDUList, the changes are written back "in place".
        assert_equal(u[0].header['rootname'], 'U2EQ0201T')
        u[0].header['rootname'] = 'abc'
        assert_equal(u[1].data[1,1], 349)
        u[1].data[1,1] = 99
        u.flush()

        # If the changes affect the size structure, e.g. adding or deleting
        # HDU(s), header was expanded or reduced beyond existing number of
        # blocks (2880 bytes in each block), or change the data size, the
        # HDUList is written to a temporary file, the original file is deleted,
        # and the temporary file is renamed to the original file name and
        # reopened in the update mode.
        # To a user, these two kinds of updating writeback seem to be the same,
        # unless the optional argument in flush or close is set to 1.
        del u[2]
        u.flush()

        # the write method in HDUList class writes the current HDUList, with
        # all changes made up to now, to a new file.  This method works the
        # same disregard the mode the HDUList was opened with.
        u.append(r[3])
        u.writeto(self.temp('test_new.fits'))

        # Remove temporary files created by this test
        u.close()


        #Another useful new HDUList method is readall.  It will "touch" the
        # data parts in all HDUs, so even if the HDUList is closed, we can
        # still operate on the data.
        r = pyfits.open(self.data('test0.fits'))
        r.readall()
        r.close()
        assert_equal(r[1].data[1,1], 315)

        # create an HDU with data only
        data = np.ones((3,5), dtype=np.float32)
        hdu = pyfits.ImageHDU(data=data, name='SCI')
        assert_equal(hdu.data.all(),
                     np.array([[ 1.,  1.,  1.,  1.,  1.],
                               [ 1.,  1.,  1.,  1.,  1.],
                               [ 1.,  1.,  1.,  1.,  1.]],
                               dtype=np.float32).all())


        # create an HDU with header and data
        # notice that the header has the right NAXIS's since it is constructed
        # with ImageHDU
        hdu2 = pyfits.ImageHDU(header=r[1].header, data=np.array([1,2],
                               dtype='int32'))

        assert_equal(str(hdu2.header.ascard[1:5]),
            "BITPIX  =                   32 / array data type                                \n"
            "NAXIS   =                    1 / number of array dimensions                     \n"
            "NAXIS1  =                    2                                                  \n"
            "PCOUNT  =                    0 / number of parameters                           ")

    def test_memory_mapping(self):
        # memory mapping
        f1 = pyfits.open(self.data('test0.fits'), memmap=1)
        f1.close()

    def test_verification_on_output(self):
        # verification on output
        # make a defect HDUList first
        with CaptureStdout() as f:
            x = pyfits.ImageHDU()
            hdu = pyfits.HDUList(x) # HDUList can take a list or one single HDU
            hdu.verify()
            assert_true(
                "HDUList's 0th element is not a primary HDU." in f.getvalue())

        with CaptureStdout() as f:
            hdu.writeto(self.temp('test_new2.fits'), 'fix')
            assert_true(
                "HDUList's 0th element is not a primary HDU.  "
                "Fixed by inserting one as 0th HDU." in f.getvalue())

    def test_section(self):
        # section testing
        fs = pyfits.open(self.data('arange.fits'))
        assert_equal(fs[0].section[3,2,5], np.array([357]))
        assert_equal(fs[0].section[3,2,:].all(),
                     np.array([352, 353, 354, 355, 356, 357, 358, 359, 360,
                               361, 362]).all())
        assert_equal(fs[0].section[3,2,4:].all(),
                     np.array([356, 357, 358, 359, 360, 361, 362]).all())
        assert_equal(fs[0].section[3,2,:8].all(),
                     np.array([352, 353, 354, 355, 356, 357, 358, 359]).all())
        assert_equal(fs[0].section[3,2,-8:8].all(),
                     np.array([355, 356, 357, 358, 359]).all())
        assert_equal(fs[0].section[3,2:5,:].all(),
                     np.array([[352, 353, 354, 355, 356, 357, 358, 359,
                                360, 361, 362],
                               [363, 364, 365, 366, 367, 368, 369, 370,
                                371, 372, 373],
                               [374, 375, 376, 377, 378, 379, 380, 381,
                                382, 383, 384]]).all())

        assert_equal(fs[0].section[3,:,:][:3,:3].all(),
                     np.array([[330, 331, 332],
                               [341, 342, 343],
                               [352, 353, 354]]).all())

        dat = fs[0].data
        assert_equal(fs[0].section[3,2:5,:8].all(), dat[3,2:5,:8].all())
        assert_equal(fs[0].section[3,2:5,3].all(), dat[3,2:5,3].all())

        assert_equal(fs[0].section[3:6,:,:][:3,:3,:3].all(),
                     np.array([[[330, 331, 332],
                                [341, 342, 343],
                                [352, 353, 354]],
                               [[440, 441, 442],
                                [451, 452, 453],
                                [462, 463, 464]],
                               [[550, 551, 552],
                                [561, 562, 563],
                                [572, 573, 574]]]).all())

        assert_equal(fs[0].section[:,:,:][:3,:2,:2].all(),
                     np.array([[[  0,   1],
                                [ 11,  12]],
                               [[110, 111],
                                [121, 122]],
                               [[220, 221],
                                [231, 232]]]).all())

        assert_equal(fs[0].section[:,2,:].all(), dat[:,2,:].all())
        assert_equal(fs[0].section[:,2:5,:].all(), dat[:,2:5,:].all())
        assert_equal(fs[0].section[3:6,3,:].all(), dat[3:6,3,:].all())
        assert_equal(fs[0].section[3:6,3:7,:].all(), dat[3:6,3:7,:].all())

    def test_section_data_square(self):
        a = np.arange(4).reshape((2, 2))
        hdu = pyfits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = pyfits.open(self.temp('test_new.fits'))
        d = hdul[0]
        dat = hdul[0].data
        assert_equal(d.section[:,:].all(), dat[:,:].all())
        assert_equal(d.section[0,:].all(), dat[0,:].all())
        assert_equal(d.section[1,:].all(), dat[1,:].all())
        assert_equal(d.section[:,0].all(), dat[:,0].all())
        assert_equal(d.section[:,1].all(), dat[:,1].all())
        assert_equal(d.section[0,0].all(), dat[0,0].all())
        assert_equal(d.section[0,1].all(), dat[0,1].all())
        assert_equal(d.section[1,0].all(), dat[1,0].all())
        assert_equal(d.section[1,1].all(), dat[1,1].all())
        assert_equal(d.section[0:1,0:1].all(), dat[0:1,0:1].all())
        assert_equal(d.section[0:2,0:1].all(), dat[0:2,0:1].all())
        assert_equal(d.section[0:1,0:2].all(), dat[0:1,0:2].all())
        assert_equal(d.section[0:2,0:2].all(), dat[0:2,0:2].all())

    def test_section_data_cube(self):
        a=np.arange(18).reshape((2,3,3))
        hdu = pyfits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul=pyfits.open(self.temp('test_new.fits'))
        d = hdul[0]
        dat = hdul[0].data
        assert_equal(d.section[:,:,:].all(), dat[:,:,:].all())
        assert_equal(d.section[:,:].all(), dat[:,:].all())
        assert_equal(d.section[:].all(), dat[:].all())
        assert_equal(d.section[0,:,:].all(), dat[0,:,:].all())
        assert_equal(d.section[1,:,:].all(), dat[1,:,:].all())
        assert_equal(d.section[0,0,:].all(), dat[0,0,:].all())
        assert_equal(d.section[0,1,:].all(), dat[0,1,:].all())
        assert_equal(d.section[0,2,:].all(), dat[0,2,:].all())
        assert_equal(d.section[1,0,:].all(), dat[1,0,:].all())
        assert_equal(d.section[1,1,:].all(), dat[1,1,:].all())
        assert_equal(d.section[1,2,:].all(), dat[1,2,:].all())
        assert_equal(d.section[0,0,0].all(), dat[0,0,0].all())
        assert_equal(d.section[0,0,1].all(), dat[0,0,1].all())
        assert_equal(d.section[0,0,2].all(), dat[0,0,2].all())
        assert_equal(d.section[0,1,0].all(), dat[0,1,0].all())
        assert_equal(d.section[0,1,1].all(), dat[0,1,1].all())
        assert_equal(d.section[0,1,2].all(), dat[0,1,2].all())
        assert_equal(d.section[0,2,0].all(), dat[0,2,0].all())
        assert_equal(d.section[0,2,1].all(), dat[0,2,1].all())
        assert_equal(d.section[0,2,2].all(), dat[0,2,2].all())
        assert_equal(d.section[1,0,0].all(), dat[1,0,0].all())
        assert_equal(d.section[1,0,1].all(), dat[1,0,1].all())
        assert_equal(d.section[1,0,2].all(), dat[1,0,2].all())
        assert_equal(d.section[1,1,0].all(), dat[1,1,0].all())
        assert_equal(d.section[1,1,1].all(), dat[1,1,1].all())
        assert_equal(d.section[1,1,2].all(), dat[1,1,2].all())
        assert_equal(d.section[1,2,0].all(), dat[1,2,0].all())
        assert_equal(d.section[1,2,1].all(), dat[1,2,1].all())
        assert_equal(d.section[1,2,2].all(), dat[1,2,2].all())
        assert_equal(d.section[:,0,0].all(), dat[:,0,0].all())
        assert_equal(d.section[:,0,1].all(), dat[:,0,1].all())
        assert_equal(d.section[:,0,2].all(), dat[:,0,2].all())
        assert_equal(d.section[:,1,0].all(), dat[:,1,0].all())
        assert_equal(d.section[:,1,1].all(), dat[:,1,1].all())
        assert_equal(d.section[:,1,2].all(), dat[:,1,2].all())
        assert_equal(d.section[:,2,0].all(), dat[:,2,0].all())
        assert_equal(d.section[:,2,1].all(), dat[:,2,1].all())
        assert_equal(d.section[:,2,2].all(), dat[:,2,2].all())
        assert_equal(d.section[0,:,0].all(), dat[0,:,0].all())
        assert_equal(d.section[0,:,1].all(), dat[0,:,1].all())
        assert_equal(d.section[0,:,2].all(), dat[0,:,2].all())
        assert_equal(d.section[1,:,0].all(), dat[1,:,0].all())
        assert_equal(d.section[1,:,1].all(), dat[1,:,1].all())
        assert_equal(d.section[1,:,2].all(), dat[1,:,2].all())
        assert_equal(d.section[:,:,0].all(), dat[:,:,0].all())
        assert_equal(d.section[:,:,1].all(), dat[:,:,1].all())
        assert_equal(d.section[:,:,2].all(), dat[:,:,2].all())
        assert_equal(d.section[:,0,:].all(), dat[:,0,:].all())
        assert_equal(d.section[:,1,:].all(), dat[:,1,:].all())
        assert_equal(d.section[:,2,:].all(), dat[:,2,:].all())

        assert_equal(d.section[:,:,0:1].all(), dat[:,:,0:1].all())
        assert_equal(d.section[:,:,0:2].all(), dat[:,:,0:2].all())
        assert_equal(d.section[:,:,0:3].all(), dat[:,:,0:3].all())
        assert_equal(d.section[:,:,1:2].all(), dat[:,:,1:2].all())
        assert_equal(d.section[:,:,1:3].all(), dat[:,:,1:3].all())
        assert_equal(d.section[:,:,2:3].all(), dat[:,:,2:3].all())
        assert_equal(d.section[0:1,0:1,0:1].all(), dat[0:1,0:1,0:1].all())
        assert_equal(d.section[0:1,0:1,0:2].all(), dat[0:1,0:1,0:2].all())
        assert_equal(d.section[0:1,0:1,0:3].all(), dat[0:1,0:1,0:3].all())
        assert_equal(d.section[0:1,0:1,1:2].all(), dat[0:1,0:1,1:2].all())
        assert_equal(d.section[0:1,0:1,1:3].all(), dat[0:1,0:1,1:3].all())
        assert_equal(d.section[0:1,0:1,2:3].all(), dat[0:1,0:1,2:3].all())
        assert_equal(d.section[0:1,0:2,0:1].all(), dat[0:1,0:2,0:1].all())
        assert_equal(d.section[0:1,0:2,0:2].all(), dat[0:1,0:2,0:2].all())
        assert_equal(d.section[0:1,0:2,0:3].all(), dat[0:1,0:2,0:3].all())
        assert_equal(d.section[0:1,0:2,1:2].all(), dat[0:1,0:2,1:2].all())
        assert_equal(d.section[0:1,0:2,1:3].all(), dat[0:1,0:2,1:3].all())
        assert_equal(d.section[0:1,0:2,2:3].all(), dat[0:1,0:2,2:3].all())
        assert_equal(d.section[0:1,0:3,0:1].all(), dat[0:1,0:3,0:1].all())
        assert_equal(d.section[0:1,0:3,0:2].all(), dat[0:1,0:3,0:2].all())
        assert_equal(d.section[0:1,0:3,0:3].all(), dat[0:1,0:3,0:3].all())
        assert_equal(d.section[0:1,0:3,1:2].all(), dat[0:1,0:3,1:2].all())
        assert_equal(d.section[0:1,0:3,1:3].all(), dat[0:1,0:3,1:3].all())
        assert_equal(d.section[0:1,0:3,2:3].all(), dat[0:1,0:3,2:3].all())
        assert_equal(d.section[0:1,1:2,0:1].all(), dat[0:1,1:2,0:1].all())
        assert_equal(d.section[0:1,1:2,0:2].all(), dat[0:1,1:2,0:2].all())
        assert_equal(d.section[0:1,1:2,0:3].all(), dat[0:1,1:2,0:3].all())
        assert_equal(d.section[0:1,1:2,1:2].all(), dat[0:1,1:2,1:2].all())
        assert_equal(d.section[0:1,1:2,1:3].all(), dat[0:1,1:2,1:3].all())
        assert_equal(d.section[0:1,1:2,2:3].all(), dat[0:1,1:2,2:3].all())
        assert_equal(d.section[0:1,1:3,0:1].all(), dat[0:1,1:3,0:1].all())
        assert_equal(d.section[0:1,1:3,0:2].all(), dat[0:1,1:3,0:2].all())
        assert_equal(d.section[0:1,1:3,0:3].all(), dat[0:1,1:3,0:3].all())
        assert_equal(d.section[0:1,1:3,1:2].all(), dat[0:1,1:3,1:2].all())
        assert_equal(d.section[0:1,1:3,1:3].all(), dat[0:1,1:3,1:3].all())
        assert_equal(d.section[0:1,1:3,2:3].all(), dat[0:1,1:3,2:3].all())
        assert_equal(d.section[1:2,0:1,0:1].all(), dat[1:2,0:1,0:1].all())
        assert_equal(d.section[1:2,0:1,0:2].all(), dat[1:2,0:1,0:2].all())
        assert_equal(d.section[1:2,0:1,0:3].all(), dat[1:2,0:1,0:3].all())
        assert_equal(d.section[1:2,0:1,1:2].all(), dat[1:2,0:1,1:2].all())
        assert_equal(d.section[1:2,0:1,1:3].all(), dat[1:2,0:1,1:3].all())
        assert_equal(d.section[1:2,0:1,2:3].all(), dat[1:2,0:1,2:3].all())
        assert_equal(d.section[1:2,0:2,0:1].all(), dat[1:2,0:2,0:1].all())
        assert_equal(d.section[1:2,0:2,0:2].all(), dat[1:2,0:2,0:2].all())
        assert_equal(d.section[1:2,0:2,0:3].all(), dat[1:2,0:2,0:3].all())
        assert_equal(d.section[1:2,0:2,1:2].all(), dat[1:2,0:2,1:2].all())
        assert_equal(d.section[1:2,0:2,1:3].all(), dat[1:2,0:2,1:3].all())
        assert_equal(d.section[1:2,0:2,2:3].all(), dat[1:2,0:2,2:3].all())
        assert_equal(d.section[1:2,0:3,0:1].all(), dat[1:2,0:3,0:1].all())
        assert_equal(d.section[1:2,0:3,0:2].all(), dat[1:2,0:3,0:2].all())
        assert_equal(d.section[1:2,0:3,0:3].all(), dat[1:2,0:3,0:3].all())
        assert_equal(d.section[1:2,0:3,1:2].all(), dat[1:2,0:3,1:2].all())
        assert_equal(d.section[1:2,0:3,1:3].all(), dat[1:2,0:3,1:3].all())
        assert_equal(d.section[1:2,0:3,2:3].all(), dat[1:2,0:3,2:3].all())
        assert_equal(d.section[1:2,1:2,0:1].all(), dat[1:2,1:2,0:1].all())
        assert_equal(d.section[1:2,1:2,0:2].all(), dat[1:2,1:2,0:2].all())
        assert_equal(d.section[1:2,1:2,0:3].all(), dat[1:2,1:2,0:3].all())
        assert_equal(d.section[1:2,1:2,1:2].all(), dat[1:2,1:2,1:2].all())
        assert_equal(d.section[1:2,1:2,1:3].all(), dat[1:2,1:2,1:3].all())
        assert_equal(d.section[1:2,1:2,2:3].all(), dat[1:2,1:2,2:3].all())
        assert_equal(d.section[1:2,1:3,0:1].all(), dat[1:2,1:3,0:1].all())
        assert_equal(d.section[1:2,1:3,0:2].all(), dat[1:2,1:3,0:2].all())
        assert_equal(d.section[1:2,1:3,0:3].all(), dat[1:2,1:3,0:3].all())
        assert_equal(d.section[1:2,1:3,1:2].all(), dat[1:2,1:3,1:2].all())
        assert_equal(d.section[1:2,1:3,1:3].all(), dat[1:2,1:3,1:3].all())
        assert_equal(d.section[1:2,1:3,2:3].all(), dat[1:2,1:3,2:3].all())

    def test_section_data_four(self):
        a = np.arange(256).reshape((4, 4, 4, 4))
        hdu = pyfits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul=pyfits.open(self.temp('test_new.fits'))
        d=hdul[0]
        dat = hdul[0].data
        assert_equal(d.section[:,:,:,:].all(), dat[:,:,:,:].all())
        assert_equal(d.section[:,:,:].all(), dat[:,:,:].all())
        assert_equal(d.section[:,:].all(), dat[:,:].all())
        assert_equal(d.section[:].all(), dat[:].all())
        assert_equal(d.section[0,:,:,:].all(), dat[0,:,:,:].all())
        assert_equal(d.section[0,:,0,:].all(), dat[0,:,0,:].all())
        assert_equal(d.section[:,:,0,:].all(), dat[:,:,0,:].all())
        assert_equal(d.section[:,1,0,:].all(), dat[:,1,0,:].all())
        assert_equal(d.section[:,:,:,1].all(), dat[:,:,:,1].all())

    def _test_comp_image(self, data, compression_type, quantize_level):
        primary_hdu = pyfits.PrimaryHDU()
        ofd = pyfits.HDUList(primary_hdu)
        chdu = pyfits.CompImageHDU(data, name='SCI',
                                   compressionType=compression_type,
                                   quantizeLevel=quantize_level)
        ofd.append(chdu)
        ofd.writeto(self.temp('test_new.fits'))
        ofd.close()
        fd = pyfits.open(self.temp('test_new.fits'))
        assert_equal(fd[1].data.all(), data.all())
        assert_equal(fd[1].header['NAXIS'], chdu.header['NAXIS'])
        assert_equal(fd[1].header['NAXIS1'], chdu.header['NAXIS1'])
        assert_equal(fd[1].header['NAXIS2'], chdu.header['NAXIS2'])
        assert_equal(fd[1].header['BITPIX'], chdu.header['BITPIX'])
        fd.close()

    def test_comp_image_rice_1(self):
        """Tests image compression with the RICE_1 algorithm."""

        self._test_comp_image(np.zeros((2, 10, 10), dtype=np.float32),
                              'RICE_1', 16)

    def test_comp_image_gzip_1(self):
        """Tests image compression with the GZIP_1 algorithm."""

        self._test_comp_image(np.zeros((2, 10, 10), dtype=np.float32),
                              'GZIP_1', -0.01)

    def test_comp_image_hcompression_1(self):
        """Tests image compression with the HCOMPRESS_1 algorithm.

        This is not a comprehensive test--just a simple round-trip test to make
        sure the code for handling HCOMPRESS_1 at least gets exercised.
        """

        assert_raises(ValueError, self._test_comp_image,
                      np.zeros((2, 10, 10), dtype=np.float32), 'HCOMPRESS_1',
                      16)
        self._test_comp_image(np.zeros((100, 100)) + 1, 'HCOMPRESS_1', 16)

    def test_disable_image_compression(self):
        with catch_warnings():
            # No warnings should be displayed in this case
            warnings.simplefilter('error')
            hdul = pyfits.open(self.data('comp.fits'),
                               disable_image_compression=True)
            # The compressed image HDU should show up as a BinTableHDU, but
            # *not* a CompImageHDU
            assert_true(isinstance(hdul[1], pyfits.BinTableHDU))
            assert_false(isinstance(hdul[1], pyfits.CompImageHDU))

    def test_do_not_scale_image_data(self):
        hdul = pyfits.open(self.data('scale.fits'),
                           do_not_scale_image_data=True)
        assert_equal(hdul[0].data.dtype, np.dtype('>i2'))
        hdul = pyfits.open(self.data('scale.fits'))
        assert_equal(hdul[0].data.dtype, np.dtype('float32'))

    def test_append_uint_data(self):
        """Test for ticket #56 (BZERO and BSCALE added in the wrong location
        when appending scaled data)
        """

        pyfits.writeto(self.temp('test_new.fits'), data=np.array([],
                       dtype='uint8'))
        d = np.zeros([100, 100]).astype('uint16')
        pyfits.append(self.temp('test_new.fits'), data=d)
        f = pyfits.open(self.temp('test_new.fits'), uint=True)
        assert_equal(f[1].data.dtype, 'uint16')

    def test_blanks(self):
        """Test image data with blank spots in it (which should show up as
        NaNs in the data array.
        """

        arr = np.zeros((10, 10), dtype=np.int32)
        # One row will be blanks
        arr[1] = 999
        hdu = pyfits.ImageHDU(data=arr)
        hdu.header.update('BLANK', 999)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = pyfits.open(self.temp('test_new.fits'))
        assert_true(np.isnan(hdul[1].data[1]).all())

    def test_bzero_with_floats(self):
        """Test use of the BZERO keyword in an image HDU containing float
        data.
        """

        arr = np.zeros((10, 10)) - 1
        hdu = pyfits.ImageHDU(data=arr)
        hdu.header.update('BZERO', 1.0)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = pyfits.open(self.temp('test_new.fits'))
        arr += 1
        assert_true((hdul[1].data == arr).all())

    def test_rewriting_large_scaled_image(self):
        """Regression test for #84"""

        hdul = pyfits.open(self.data('fixed-1890.fits'))
        orig_data = hdul[0].data
        hdul.writeto(self.temp('test_new.fits'), clobber=True)
        hdul.close()
        hdul = pyfits.open(self.temp('test_new.fits'))
        assert_true((hdul[0].data == orig_data).all())
        hdul.close()

        # Just as before, but this time don't touch hdul[0].data before writing
        # back out--this is the case that failed in #84
        hdul = pyfits.open(self.data('fixed-1890.fits'))
        hdul.writeto(self.temp('test_new.fits'), clobber=True)
        hdul.close()
        hdul = pyfits.open(self.temp('test_new.fits'))
        assert_true((hdul[0].data == orig_data).all())
        hdul.close()
