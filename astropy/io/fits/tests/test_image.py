import os
import warnings

import numpy as np

from ....io import fits
from ....tests.helper import pytest, raises

from . import FitsTestCase


class TestImageFunctions(FitsTestCase):
    def test_card_constructor_default_args(self):
        """Test the constructor with default argument values."""

        c = fits.Card()
        assert '' == c.key

    def test_constructor_name_arg(self):
        """Like the test of the same name in test_table.py"""

        hdu = fits.ImageHDU()
        assert hdu.name == ''
        assert 'EXTNAME' not in hdu.header
        hdu.name = 'FOO'
        assert hdu.name == 'FOO'
        assert hdu.header['EXTNAME'] == 'FOO'

        # Passing name to constructor
        hdu = fits.ImageHDU(name='FOO')
        assert hdu.name == 'FOO'
        assert hdu.header['EXTNAME'] == 'FOO'

        # And overriding a header with a different extname
        hdr = fits.Header()
        hdr.update('EXTNAME', 'EVENTS')
        hdu = fits.ImageHDU(header=hdr, name='FOO')
        assert hdu.name == 'FOO'
        assert hdu.header['EXTNAME'] == 'FOO'

    def test_fromstring_set_attribute_ascardimage(self):
        """Test fromstring() which will return a new card."""

        c = fits.Card('abc', 99).fromstring('xyz     = 100')
        assert 100 == c.value

        # test set attribute and  ascardimage() using the most updated attributes
        c.value = 200
        assert (c.ascardimage() ==
                "XYZ     =                  200                                                  ")

    def test_string(self):
        """Test string value"""

        c = fits.Card('abc', '<8 ch')
        assert (str(c) ==
                "ABC     = '<8 ch   '                                                            ")
        c = fits.Card('nullstr', '')
        assert(str(c) ==
               "NULLSTR = ''                                                                    ")

    def test_boolean_value_card(self):
        """Boolean value card"""

        c = fits.Card("abc", True)
        assert (str(c) ==
                "ABC     =                    T                                                  ")

        c = fits.Card.fromstring('abc     = F')
        assert c.value == False

    def test_long_integer_number(self):
        """long integer number"""

        c = fits.Card('long_int', -467374636747637647347374734737437)
        assert (str(c) ==
                "LONG_INT= -467374636747637647347374734737437                                    ")

    def test_floating_point_number(self):
        """ floating point number"""

        c = fits.Card('floatnum', -467374636747637647347374734737437.)

        if (str(c) != "FLOATNUM= -4.6737463674763E+32                                                  " and
            str(c) != "FLOATNUM= -4.6737463674763E+032                                                 "):
            assert (str(c) ==
                    "FLOATNUM= -4.6737463674763E+32                                                  ")

    def test_complex_value(self):
        """complex value"""

        c = fits.Card('abc',
                        1.2345377437887837487e88+6324767364763746367e-33j)

        if (str(c) != "ABC     = (1.23453774378878E+88, 6.32476736476374E-15)                          " and
            str(c) != "ABC     = (1.2345377437887E+088, 6.3247673647637E-015)                          "):
            assert (str(c),
                    "ABC     = (1.23453774378878E+88, 6.32476736476374E-15)                          ")

    def test_card_image_constructed_too_long(self):
        # card image constructed from key/value/comment is too long (non-string
        # value)
        c = fits.Card('abc', 9, 'abcde'*20)
        assert (str(c) ==
                "ABC     =                    9 / abcdeabcdeabcdeabcdeabcdeabcdeabcdeabcdeabcdeab")
        c = fits.Card('abc', 'a'*68, 'abcdefg')
        assert (str(c) ==
                "ABC     = 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'")

    def test_constructor_filter_illegal_data_structures(self):
        # the constrctor will filter out illegal data structures...
        pytest.raises(ValueError, fits.Card, ('abc',), {'value': (2, 3)})
        pytest.raises(ValueError, fits.Card, 'key', [], 'comment')

    @raises(ValueError)
    def test_keyword_too_long(self):
        #keywords too long
        fits.Card('abcdefghi', 'long')

    @raises(ValueError)
    def test_illegal_characters_in_key(self):
        # will not allow illegal characters in key when using constructor
        fits.Card('abc+', 9)

    # TODO: What sort of exception should this raise?
    @raises(Exception)
    def test_ascardimage_verifies_the_comment_string_to_be_ascii_text(self):
        # the ascardimage() verifies the comment string to be ASCII text
        c = fits.Card.fromstring('abc     = +  2.1   e + 12 / abcde\0')
        c.ascardimage()

    def test_commentary_cards(self):
        # commentary cards
        c = fits.Card("history",
                        "A commentary card's value has no quotes around it.")
        assert (str(c) ==
                "HISTORY A commentary card's value has no quotes around it.                      ")
        c = fits.Card("comment",
                        "A commentary card has no comment.", "comment")
        assert (str(c) ==
                "COMMENT A commentary card has no comment.                                       ")

    def test_commentary_card_created_by_fromstring(self):
        # commentary card created by fromstring()
        c = fits.Card.fromstring("COMMENT card has no comments. / text after slash is still part of the value.")
        assert (c.value ==
                'card has no comments. / text after slash is still part of the value.')
        assert c.comment == ''

    def test_commentary_card_will_not_parse_numerical_value(self):
        # commentary card will not parse the numerical value
        c = fits.Card.fromstring("history  (1, 2)")
        assert (str(c.ascardimage()) ==
                "HISTORY  (1, 2)                                                                 ")

    def test_equal_sign_after_column8(self):
        # equal sign after column 8 of a commentary card will be part ofthe string value
        c = fits.Card.fromstring("history =   (1, 2)")
        assert (str(c.ascardimage()) ==
                "HISTORY =   (1, 2)                                                              ")

    def test_specify_undefined_value(self):
        # this is how to specify an undefined value
        c = fits.Card("undef", fits.card.UNDEFINED)
        assert (str(c) ==
                     "UNDEF   =                                                                       ")

    def test_complex_number_using_string_input(self):
        # complex number using string input
        c = fits.Card.fromstring('abc     = (8, 9)')
        assert (str(c.ascardimage()) ==
                "ABC     =               (8, 9)                                                  ")

    def test_fixable_non_standard_fits_card(self):
        # fixable non-standard FITS card will keep the original format
        c = fits.Card.fromstring('abc     = +  2.1   e + 12')
        assert c.value == 2100000000000.0
        assert (str(c.ascardimage()) ==
                "ABC     =             +2.1E+12                                                  ")

    def test_fixable_non_fsc(self):
        # fixable non-FSC: if the card is not parsable, it's value will be
        # assumed
        # to be a string and everything after the first slash will be comment
        c = fits.Card.fromstring("no_quote=  this card's value has no quotes / let's also try the comment")
        assert (str(c.ascardimage()) ==
                "NO_QUOTE= 'this card''s value has no quotes' / let's also try the comment       ")

    def test_undefined_value_using_string_input(self):
        # undefined value using string input
        c = fits.Card.fromstring('abc     =    ')
        assert (str(c.ascardimage()) ==
                "ABC     =                                                                       ")

    def test_misalocated_equal_sign(self):
        # test mislocated "=" sign
        c = fits.Card.fromstring('xyz= 100')
        assert c.key == 'XYZ'
        assert c.value == 100
        assert (str(c.ascardimage()) ==
                "XYZ     =                  100                                                  ")

    def test_equal_only_up_to_column_10(self):
        # the test of "=" location is only up to column 10
        c = fits.Card.fromstring("histo       =   (1, 2)")
        assert (str(c.ascardimage()) ==
                "HISTO   = '=   (1, 2)'                                                          ")
        c = fits.Card.fromstring("   history          (1, 2)")
        assert (str(c.ascardimage()) ==
                "HISTO   = 'ry          (1, 2)'                                                  ")

    def test_verification(self, capsys):
        # verification
        c = fits.Card.fromstring('abc= a6')
        c.verify()
        out, err = capsys.readouterr()
        assert ('Card image is not FITS standard (equal sign not at column 8).'
                in err)
        assert (str(c) ==
                "abc= a6                                                                         ")

    def test_fix(self, capsys):
        c = fits.Card.fromstring('abc= a6')
        c.verify('fix')
        out, err = capsys.readouterr()
        assert 'Fixed card to be FITS standard.: ABC' in err
        assert (str(c) ==
                "ABC     = 'a6      '                                                            ")

    def test_long_string_value(self):
        # test long string value
        c = fits.Card('abc', 'long string value '*10, 'long comment '*10)
        assert (str(c) ==
                "ABC     = 'long string value long string value long string value long string &' "
                "CONTINUE  'value long string value long string value long string value long &'  "
                "CONTINUE  'string value long string value long string value &'                  "
                "CONTINUE  '&' / long comment long comment long comment long comment long        "
                "CONTINUE  '&' / comment long comment long comment long comment long comment     "
                "CONTINUE  '&' / long comment                                                    ")

    def test_long_string_from_file(self):
        c = fits.Card('abc', 'long string value '*10, 'long comment '*10)
        hdu = fits.PrimaryHDU()
        hdu.header.ascard.append(c)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        c = hdul[0].header.ascard['abc']
        hdul.close()
        assert (str(c),
                "ABC     = 'long string value long string value long string value long string &' "
                "CONTINUE  'value long string value long string value long string value long &'  "
                "CONTINUE  'string value long string value long string value &'                  "
                "CONTINUE  '&' / long comment long comment long comment long comment long        "
                "CONTINUE  '&' / comment long comment long comment long comment long comment     "
                "CONTINUE  '&' / long comment                                                    ")


    def test_word_in_long_string_too_long(self):
        # if a word in a long string is too long, it will be cut in the middle
        c = fits.Card('abc', 'longstringvalue'*10, 'longcomment'*10)
        assert (str(c) ==
                "ABC     = 'longstringvaluelongstringvaluelongstringvaluelongstringvaluelongstr&'"
                "CONTINUE  'ingvaluelongstringvaluelongstringvaluelongstringvaluelongstringvalu&'"
                "CONTINUE  'elongstringvalue&'                                                   "
                "CONTINUE  '&' / longcommentlongcommentlongcommentlongcommentlongcommentlongcomme"
                "CONTINUE  '&' / ntlongcommentlongcommentlongcommentlongcomment                  ")

    def test_long_string_value_via_fromstring(self):
        # long string value via fromstring() method
        c = fits.Card.fromstring(
            fits.card._pad("abc     = 'longstring''s testing  &  ' / comments in line 1") +
            fits.card._pad("continue  'continue with long string but without the ampersand at the end' /") +
            fits.card._pad("continue  'continue must have string value (with quotes)' / comments with ''. "))
        assert (str(c.ascardimage()) ==
                "ABC     = 'longstring''s testing  continue with long string but without the &'  "
                "CONTINUE  'ampersand at the endcontinue must have string value (with quotes)&'  "
                "CONTINUE  '&' / comments in line 1 comments with ''.                            ")

    def test_hierarch_card(self):
        c = fits.Card('hierarch abcdefghi', 10)
        assert (str(c.ascardimage()) ==
                "HIERARCH abcdefghi = 10                                                         ")
        c = fits.Card('HIERARCH ESO INS SLIT2 Y1FRML', 'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)')
        assert (str(c.ascardimage()) ==
                "HIERARCH ESO INS SLIT2 Y1FRML= 'ENC=OFFSET+RESOL*acos((WID-(MAX+MIN))/(MAX-MIN)'")


    # TODO: What sort of exception should this raise?
    @raises(Exception)
    def test_open(self):
        # The function "open" reads a FITS file into an HDUList object.  There
        # are three modes to open: "readonly" (the default), "append", and
        # "update".

        # Open a file read-only (the default mode), the content of the FITS
        # file are read into memory.
        r = fits.open(self.data('test0.fits')) # readonly

        # data parts are latent instantiation, so if we close the HDUList
        # without touching data, data can not be accessed.
        r.close()
        r[1].data[:2,:2]

    def test_open_2(self):
        r = fits.open(self.data('test0.fits'))

        info = [(0, 'PRIMARY', 'PrimaryHDU', 138, (), 'int16', '')] + \
               [(x, 'SCI', 'ImageHDU', 61, (40, 40), 'int16', '')
                for x in range(1, 5)]

        try:
            assert r.info(output=False) == info
        finally:
            r.close()

    def test_io_manipulation(self):
        # Get a keyword value.  An extension can be referred by name or by
        # number.  Both extension and keyword names are case insensitive.
        r = fits.open(self.data('test0.fits'))
        assert r['primary'].header['naxis'] == 0
        assert r[0].header['naxis'] == 0

        # If there are more than one extension with the same EXTNAME value, the
        # EXTVER can be used (as the second argument) to distinguish the
        # extension.
        assert r['sci',1].header['detector'] == 1

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
            assert (str(r[0].header.ascard[-3:]) ==
                    "EXPFLAG = 'NORMAL            ' / Exposure interruption indicator                \n"
                    "FILENAME= 'vtest3.fits'        / File name                                      \n"
                    "XXX     =            1.234E+56                                                  ")

        # rename a keyword
        r[0].header.rename_key('filename', 'fname')
        pytest.raises(ValueError, r[0].header.rename_key, 'fname', 'history')

        pytest.raises(ValueError, r[0].header.rename_key, 'fname', 'simple')
        r[0].header.rename_key('fname', 'filename')

        # get a subsection of data
        assert (r[2].data[:3,:3] ==
                np.array([[349, 349, 348],
                          [349, 349, 347],
                          [347, 350, 349]], dtype=np.int16)).all()

        # We can create a new FITS file by opening a new file with "append"
        # mode.
        n=fits.open(self.temp('test_new.fits'), mode='append')

        # Append the primary header and the 2nd extension to the new file.
        n.append(r[0])
        n.append(r[2])

        # The flush method will write the current HDUList object back to the
        # newly created file on disk.  The HDUList is still open and can be
        # further operated.
        n.flush()
        assert n[1].data[1,1] == 349

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
        a = fits.open(self.temp('test_append.fits'), mode='append')

        # The above change did not take effect since this was made after the
        # flush().
        assert a[1].data[1,1] == 349

        a.append(r[1])
        a.close()

        # When changes are made to an HDUList which was opened with "update"
        # mode, they will be written back to the original file when a
        # flush/close is called.
        os.rename(self.temp('test_append.fits'), self.temp('test_update.fits'))

        u = fits.open(self.temp('test_update.fits'), mode='update')

        # When the changes do not alter the size structures of the original (or
        # since last flush) HDUList, the changes are written back "in place".
        assert u[0].header['rootname'] == 'U2EQ0201T'
        u[0].header['rootname'] = 'abc'
        assert u[1].data[1,1] == 349
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
        r = fits.open(self.data('test0.fits'))
        r.readall()
        r.close()
        assert r[1].data[1,1] == 315

        # create an HDU with data only
        data = np.ones((3,5), dtype=np.float32)
        hdu = fits.ImageHDU(data=data, name='SCI')
        assert (hdu.data ==
                np.array([[ 1.,  1.,  1.,  1.,  1.],
                          [ 1.,  1.,  1.,  1.,  1.],
                          [ 1.,  1.,  1.,  1.,  1.]], dtype=np.float32)).all()


        # create an HDU with header and data
        # notice that the header has the right NAXIS's since it is constructed
        # with ImageHDU
        hdu2 = fits.ImageHDU(header=r[1].header, data=np.array([1,2],
                               dtype='int32'))

        assert (str(hdu2.header.ascard[1:5]) ==
                "BITPIX  =                   32 / array data type                                \n"
                "NAXIS   =                    1 / number of array dimensions                     \n"
                "NAXIS1  =                    2                                                  \n"
               "PCOUNT  =                    0 / number of parameters                           ")

    def test_memory_mapping(self):
        # memory mapping
        f1 = fits.open(self.data('test0.fits'), memmap=1)
        f1.close()

    def test_verification_on_output(self, capsys):
        # verification on output
        # make a defect HDUList first
        x = fits.ImageHDU()
        hdu = fits.HDUList(x) # HDUList can take a list or one single HDU
        hdu.verify()
        out, err = capsys.readouterr()
        assert "HDUList's 0th element is not a primary HDU." in err

        hdu.writeto(self.temp('test_new2.fits'), 'fix')
        out, err = capsys.readouterr()
        assert ("HDUList's 0th element is not a primary HDU.  "
                "Fixed by inserting one as 0th HDU." in err)

    def test_section(self):
        # section testing
        fs = fits.open(self.data('arange.fits'))
        assert fs[0].section[3,2,5] == np.array([357])
        assert (fs[0].section[3,2,:] ==
                np.array([352, 353, 354, 355, 356, 357, 358, 359, 360, 361,
                          362])).all()
        assert (fs[0].section[3,2,4:] ==
                np.array([356, 357, 358, 359, 360, 361, 362])).all()
        assert (fs[0].section[3,2,:8] ==
                np.array([352, 353, 354, 355, 356, 357, 358, 359])).all()
        assert (fs[0].section[3,2,-8:8] ==
                np.array([355, 356, 357, 358, 359])).all()
        assert (fs[0].section[3,2:5,:] ==
                np.array([[352, 353, 354, 355, 356, 357, 358, 359, 360, 361,
                           362],
                          [363, 364, 365, 366, 367, 368, 369, 370, 371, 372,
                           373],
                          [374, 375, 376, 377, 378, 379, 380, 381, 382, 383,
                           384]])).all()

        assert (fs[0].section[3,:,:][:3,:3] ==
                np.array([[330, 331, 332],
                          [341, 342, 343],
                          [352, 353, 354]])).all()

        dat = fs[0].data
        assert (fs[0].section[3,2:5,:8] == dat[3,2:5,:8]).all()
        assert (fs[0].section[3,2:5,3] == dat[3,2:5,3]).all()

        assert (fs[0].section[3:6,:,:][:3,:3,:3] ==
                np.array([[[330, 331, 332],
                           [341, 342, 343],
                           [352, 353, 354]],
                          [[440, 441, 442],
                           [451, 452, 453],
                           [462, 463, 464]],
                          [[550, 551, 552],
                           [561, 562, 563],
                           [572, 573, 574]]])).all()

        assert (fs[0].section[:,:,:][:3,:2,:2] ==
                np.array([[[  0,   1],
                           [ 11,  12]],
                          [[110, 111],
                           [121, 122]],
                          [[220, 221],
                           [231, 232]]])).all()

        assert (fs[0].section[:,2,:] == dat[:,2,:]).all()
        assert (fs[0].section[:,2:5,:] == dat[:,2:5,:]).all()
        assert (fs[0].section[3:6,3,:] == dat[3:6,3,:]).all()
        assert (fs[0].section[3:6,3:7,:] == dat[3:6,3:7,:]).all()

    def test_section_data_square(self):
        a = np.arange(4).reshape((2, 2))
        hdu = fits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        d = hdul[0]
        dat = hdul[0].data
        assert (d.section[:,:] == dat[:,:]).all()
        assert (d.section[0,:] == dat[0,:]).all()
        assert (d.section[1,:] == dat[1,:]).all()
        assert (d.section[:,0] == dat[:,0]).all()
        assert (d.section[:,1] == dat[:,1]).all()
        assert (d.section[0,0] == dat[0,0]).all()
        assert (d.section[0,1] == dat[0,1]).all()
        assert (d.section[1,0] == dat[1,0]).all()
        assert (d.section[1,1] == dat[1,1]).all()
        assert (d.section[0:1,0:1] == dat[0:1,0:1]).all()
        assert (d.section[0:2,0:1] == dat[0:2,0:1]).all()
        assert (d.section[0:1,0:2] == dat[0:1,0:2]).all()
        assert (d.section[0:2,0:2] == dat[0:2,0:2]).all()

    def test_section_data_cube(self):
        a=np.arange(18).reshape((2,3,3))
        hdu = fits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul=fits.open(self.temp('test_new.fits'))
        d = hdul[0]
        dat = hdul[0].data
        assert (d.section[:,:,:] == dat[:,:,:]).all()
        assert (d.section[:,:] == dat[:,:]).all()
        assert d.section[:].all() == dat[:].all()
        assert (d.section[0,:,:] == dat[0,:,:]).all()
        assert (d.section[1,:,:] == dat[1,:,:]).all()
        assert (d.section[0,0,:] == dat[0,0,:]).all()
        assert (d.section[0,1,:] == dat[0,1,:]).all()
        assert (d.section[0,2,:] == dat[0,2,:]).all()
        assert (d.section[1,0,:] == dat[1,0,:]).all()
        assert (d.section[1,1,:] == dat[1,1,:]).all()
        assert (d.section[1,2,:] == dat[1,2,:]).all()
        assert (d.section[0,0,0] == dat[0,0,0]).all()
        assert (d.section[0,0,1] == dat[0,0,1]).all()
        assert (d.section[0,0,2] == dat[0,0,2]).all()
        assert (d.section[0,1,0] == dat[0,1,0]).all()
        assert (d.section[0,1,1] == dat[0,1,1]).all()
        assert (d.section[0,1,2] == dat[0,1,2]).all()
        assert (d.section[0,2,0] == dat[0,2,0]).all()
        assert (d.section[0,2,1] == dat[0,2,1]).all()
        assert (d.section[0,2,2] == dat[0,2,2]).all()
        assert (d.section[1,0,0] == dat[1,0,0]).all()
        assert (d.section[1,0,1] == dat[1,0,1]).all()
        assert (d.section[1,0,2] == dat[1,0,2]).all()
        assert (d.section[1,1,0] == dat[1,1,0]).all()
        assert (d.section[1,1,1] == dat[1,1,1]).all()
        assert (d.section[1,1,2] == dat[1,1,2]).all()
        assert (d.section[1,2,0] == dat[1,2,0]).all()
        assert (d.section[1,2,1] == dat[1,2,1]).all()
        assert (d.section[1,2,2] == dat[1,2,2]).all()
        assert (d.section[:,0,0] == dat[:,0,0]).all()
        assert (d.section[:,0,1] == dat[:,0,1]).all()
        assert (d.section[:,0,2] == dat[:,0,2]).all()
        assert (d.section[:,1,0] == dat[:,1,0]).all()
        assert (d.section[:,1,1] == dat[:,1,1]).all()
        assert (d.section[:,1,2] == dat[:,1,2]).all()
        assert (d.section[:,2,0] == dat[:,2,0]).all()
        assert (d.section[:,2,1] == dat[:,2,1]).all()
        assert (d.section[:,2,2] == dat[:,2,2]).all()
        assert (d.section[0,:,0] == dat[0,:,0]).all()
        assert (d.section[0,:,1] == dat[0,:,1]).all()
        assert (d.section[0,:,2] == dat[0,:,2]).all()
        assert (d.section[1,:,0] == dat[1,:,0]).all()
        assert (d.section[1,:,1] == dat[1,:,1]).all()
        assert (d.section[1,:,2] == dat[1,:,2]).all()
        assert (d.section[:,:,0] == dat[:,:,0]).all()
        assert (d.section[:,:,1] == dat[:,:,1]).all()
        assert (d.section[:,:,2] == dat[:,:,2]).all()
        assert (d.section[:,0,:] == dat[:,0,:]).all()
        assert (d.section[:,1,:] == dat[:,1,:]).all()
        assert (d.section[:,2,:] == dat[:,2,:]).all()

        assert (d.section[:,:,0:1] == dat[:,:,0:1]).all()
        assert (d.section[:,:,0:2] == dat[:,:,0:2]).all()
        assert (d.section[:,:,0:3] == dat[:,:,0:3]).all()
        assert (d.section[:,:,1:2] == dat[:,:,1:2]).all()
        assert (d.section[:,:,1:3] == dat[:,:,1:3]).all()
        assert (d.section[:,:,2:3] == dat[:,:,2:3]).all()
        assert (d.section[0:1,0:1,0:1] == dat[0:1,0:1,0:1]).all()
        assert (d.section[0:1,0:1,0:2] == dat[0:1,0:1,0:2]).all()
        assert (d.section[0:1,0:1,0:3] == dat[0:1,0:1,0:3]).all()
        assert (d.section[0:1,0:1,1:2] == dat[0:1,0:1,1:2]).all()
        assert (d.section[0:1,0:1,1:3] == dat[0:1,0:1,1:3]).all()
        assert (d.section[0:1,0:1,2:3] == dat[0:1,0:1,2:3]).all()
        assert (d.section[0:1,0:2,0:1] == dat[0:1,0:2,0:1]).all()
        assert (d.section[0:1,0:2,0:2] == dat[0:1,0:2,0:2]).all()
        assert (d.section[0:1,0:2,0:3] == dat[0:1,0:2,0:3]).all()
        assert (d.section[0:1,0:2,1:2] == dat[0:1,0:2,1:2]).all()
        assert (d.section[0:1,0:2,1:3] == dat[0:1,0:2,1:3]).all()
        assert (d.section[0:1,0:2,2:3] == dat[0:1,0:2,2:3]).all()
        assert (d.section[0:1,0:3,0:1] == dat[0:1,0:3,0:1]).all()
        assert (d.section[0:1,0:3,0:2] == dat[0:1,0:3,0:2]).all()
        assert (d.section[0:1,0:3,0:3] == dat[0:1,0:3,0:3]).all()
        assert (d.section[0:1,0:3,1:2] == dat[0:1,0:3,1:2]).all()
        assert (d.section[0:1,0:3,1:3] == dat[0:1,0:3,1:3]).all()
        assert (d.section[0:1,0:3,2:3] == dat[0:1,0:3,2:3]).all()
        assert (d.section[0:1,1:2,0:1] == dat[0:1,1:2,0:1]).all()
        assert (d.section[0:1,1:2,0:2] == dat[0:1,1:2,0:2]).all()
        assert (d.section[0:1,1:2,0:3] == dat[0:1,1:2,0:3]).all()
        assert (d.section[0:1,1:2,1:2] == dat[0:1,1:2,1:2]).all()
        assert (d.section[0:1,1:2,1:3] == dat[0:1,1:2,1:3]).all()
        assert (d.section[0:1,1:2,2:3] == dat[0:1,1:2,2:3]).all()
        assert (d.section[0:1,1:3,0:1] == dat[0:1,1:3,0:1]).all()
        assert (d.section[0:1,1:3,0:2] == dat[0:1,1:3,0:2]).all()
        assert (d.section[0:1,1:3,0:3] == dat[0:1,1:3,0:3]).all()
        assert (d.section[0:1,1:3,1:2] == dat[0:1,1:3,1:2]).all()
        assert (d.section[0:1,1:3,1:3] == dat[0:1,1:3,1:3]).all()
        assert (d.section[0:1,1:3,2:3] == dat[0:1,1:3,2:3]).all()
        assert (d.section[1:2,0:1,0:1] == dat[1:2,0:1,0:1]).all()
        assert (d.section[1:2,0:1,0:2] == dat[1:2,0:1,0:2]).all()
        assert (d.section[1:2,0:1,0:3] == dat[1:2,0:1,0:3]).all()
        assert (d.section[1:2,0:1,1:2] == dat[1:2,0:1,1:2]).all()
        assert (d.section[1:2,0:1,1:3] == dat[1:2,0:1,1:3]).all()
        assert (d.section[1:2,0:1,2:3] == dat[1:2,0:1,2:3]).all()
        assert (d.section[1:2,0:2,0:1] == dat[1:2,0:2,0:1]).all()
        assert (d.section[1:2,0:2,0:2] == dat[1:2,0:2,0:2]).all()
        assert (d.section[1:2,0:2,0:3] == dat[1:2,0:2,0:3]).all()
        assert (d.section[1:2,0:2,1:2] == dat[1:2,0:2,1:2]).all()
        assert (d.section[1:2,0:2,1:3] == dat[1:2,0:2,1:3]).all()
        assert (d.section[1:2,0:2,2:3] == dat[1:2,0:2,2:3]).all()
        assert (d.section[1:2,0:3,0:1] == dat[1:2,0:3,0:1]).all()
        assert (d.section[1:2,0:3,0:2] == dat[1:2,0:3,0:2]).all()
        assert (d.section[1:2,0:3,0:3] == dat[1:2,0:3,0:3]).all()
        assert (d.section[1:2,0:3,1:2] == dat[1:2,0:3,1:2]).all()
        assert (d.section[1:2,0:3,1:3] == dat[1:2,0:3,1:3]).all()
        assert (d.section[1:2,0:3,2:3] == dat[1:2,0:3,2:3]).all()
        assert (d.section[1:2,1:2,0:1] == dat[1:2,1:2,0:1]).all()
        assert (d.section[1:2,1:2,0:2] == dat[1:2,1:2,0:2]).all()
        assert (d.section[1:2,1:2,0:3] == dat[1:2,1:2,0:3]).all()
        assert (d.section[1:2,1:2,1:2] == dat[1:2,1:2,1:2]).all()
        assert (d.section[1:2,1:2,1:3] == dat[1:2,1:2,1:3]).all()
        assert (d.section[1:2,1:2,2:3] == dat[1:2,1:2,2:3]).all()
        assert (d.section[1:2,1:3,0:1] == dat[1:2,1:3,0:1]).all()
        assert (d.section[1:2,1:3,0:2] == dat[1:2,1:3,0:2]).all()
        assert (d.section[1:2,1:3,0:3] == dat[1:2,1:3,0:3]).all()
        assert (d.section[1:2,1:3,1:2] == dat[1:2,1:3,1:2]).all()
        assert (d.section[1:2,1:3,1:3] == dat[1:2,1:3,1:3]).all()
        assert (d.section[1:2,1:3,2:3] == dat[1:2,1:3,2:3]).all()

    def test_section_data_four(self):
        a = np.arange(256).reshape((4, 4, 4, 4))
        hdu = fits.PrimaryHDU(a)
        hdu.writeto(self.temp('test_new.fits'))

        hdul=fits.open(self.temp('test_new.fits'))
        d=hdul[0]
        dat = hdul[0].data
        assert (d.section[:,:,:,:] == dat[:,:,:,:]).all()
        assert (d.section[:,:,:] == dat[:,:,:]).all()
        assert (d.section[:,:] == dat[:,:]).all()
        assert d.section[:].all() == dat[:].all()
        assert (d.section[0,:,:,:] == dat[0,:,:,:]).all()
        assert (d.section[0,:,0,:] == dat[0,:,0,:]).all()
        assert (d.section[:,:,0,:] == dat[:,:,0,:]).all()
        assert (d.section[:,1,0,:] == dat[:,1,0,:]).all()
        assert (d.section[:,:,:,1] == dat[:,:,:,1]).all()

    def _test_comp_image(self, data, compression_type, quantize_level):
        primary_hdu = fits.PrimaryHDU()
        ofd = fits.HDUList(primary_hdu)
        chdu = fits.CompImageHDU(data, name='SCI',
                                   compressionType=compression_type,
                                   quantizeLevel=quantize_level)
        ofd.append(chdu)
        ofd.writeto(self.temp('test_new.fits'))
        ofd.close()
        fd = fits.open(self.temp('test_new.fits'))
        assert fd[1].data.all() == data.all()
        assert fd[1].header['NAXIS'] == chdu.header['NAXIS']
        assert fd[1].header['NAXIS1'] == chdu.header['NAXIS1']
        assert fd[1].header['NAXIS2'] == chdu.header['NAXIS2']
        assert fd[1].header['BITPIX'] == chdu.header['BITPIX']
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

        pytest.raises(ValueError, self._test_comp_image,
                      np.zeros((2, 10, 10), dtype=np.float32), 'HCOMPRESS_1',
                      16)
        self._test_comp_image(np.zeros((100, 100)) + 1, 'HCOMPRESS_1', 16)

    def test_disable_image_compression(self):
        with warnings.catch_warnings():
            # No warnings should be displayed in this case
            warnings.simplefilter('error')
            hdul = fits.open(self.data('comp.fits'),
                               disable_image_compression=True)
            # The compressed image HDU should show up as a BinTableHDU, but
            # *not* a CompImageHDU
            assert isinstance(hdul[1], fits.BinTableHDU)
            assert not isinstance(hdul[1], fits.CompImageHDU)

    def test_do_not_scale_image_data(self):
        hdul = fits.open(self.data('scale.fits'),
                           do_not_scale_image_data=True)
        assert hdul[0].data.dtype == np.dtype('>i2')
        hdul = fits.open(self.data('scale.fits'))
        assert hdul[0].data.dtype == np.dtype('float32')

    def test_append_uint_data(self):
        """Test for ticket #56 (BZERO and BSCALE added in the wrong location
        when appending scaled data)
        """

        fits.writeto(self.temp('test_new.fits'), data=np.array([],
                       dtype='uint8'))
        d = np.zeros([100, 100]).astype('uint16')
        fits.append(self.temp('test_new.fits'), data=d)
        f = fits.open(self.temp('test_new.fits'), uint=True)
        assert f[1].data.dtype == 'uint16'

    def test_blanks(self):
        """Test image data with blank spots in it (which should show up as
        NaNs in the data array.
        """

        arr = np.zeros((10, 10), dtype=np.int32)
        # One row will be blanks
        arr[1] = 999
        hdu = fits.ImageHDU(data=arr)
        hdu.header.update('BLANK', 999)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        assert np.isnan(hdul[1].data[1]).all()

    def test_bzero_with_floats(self):
        """Test use of the BZERO keyword in an image HDU containing float
        data.
        """

        arr = np.zeros((10, 10)) - 1
        hdu = fits.ImageHDU(data=arr)
        hdu.header.update('BZERO', 1.0)
        hdu.writeto(self.temp('test_new.fits'))

        hdul = fits.open(self.temp('test_new.fits'))
        arr += 1
        assert (hdul[1].data == arr).all()

    def test_rewriting_large_scaled_image(self):
        """Regression test for #84"""

        hdul = fits.open(self.data('fixed-1890.fits'))
        orig_data = hdul[0].data
        hdul.writeto(self.temp('test_new.fits'), clobber=True)
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        assert (hdul[0].data == orig_data).all()
        hdul.close()

        # Just as before, but this time don't touch hdul[0].data before writing
        # back out--this is the case that failed in #84
        hdul = fits.open(self.data('fixed-1890.fits'))
        hdul.writeto(self.temp('test_new.fits'), clobber=True)
        hdul.close()
        hdul = fits.open(self.temp('test_new.fits'))
        assert (hdul[0].data == orig_data).all()
        hdul.close()
