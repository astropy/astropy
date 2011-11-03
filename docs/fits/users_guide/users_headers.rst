.. currentmodule:: pyfits.core

************
FITS Headers
************

In the next three chapters, more detailed information as well as examples will
be explained for manipulating the header, the image data, and the table data
respectively.


Header of an HDU
================

Every HDU normally has two components: header and data. In PyFITS these two
components are accessed through the two attributes of the HDU, ``.header`` and
``.data``.

While an HDU may have empty data, i.e. the .data attribute is None, any HDU
will always have a header. When an HDU is created with a constructor, e.g.
``hdu = PrimaryHDU(data, header)``, the user may supply the header value from
an existing HDU's header and the data value from  a numpy array. If the
defaults (``None``) are used, the new HDU will have the minimal require
keyword:

    >>> hdu = pyfits.PrimaryHDU()
    >>> print hdu.header.ascardlist() # show the keywords
    SIMPLE = T / conforms to FITS standard
    BITPIX = 8 / array data type
    NAXIS  = 0 / number of array dimensions
    EXTEND = T

A user can use any header and any data to construct a new HDU. PyFITS will
strip the required keywords from the input header first and then add back the
required keywords compatible to the new HDU. So, a user can use a table HDU's
header to construct an image HDU and vice versa. The constructor will also
ensure the data type and dimension information in the header agree with the
data.


The Header Attribute
====================


Value Access and Updating
-------------------------

As shown in the Quick Tutorial, keyword values can be accessed via keyword name
or index of an HDU's header attribute. Here is a quick summary:

    >>> hdulist = pyfits.open('input.fits') # open a FITS file
    >>> prihdr = hdulist[0].header # the primary HDU header
    >>> print prihdr[3] # get the 4th keyword's value
    10
    >>> prihdr[3] = 20 # change it's value
    >>> print prihdr['darkcorr'] # get the value of the keyword 'darkcorr'
    'OMIT'
    >>> prihdr['darkcorr'] = 'PERFORM' # change darkcorr's value

When reference by the keyword name, it is case insensitive. Thus,
prihdr['abc'], prihdr['ABC'], or prihdr['aBc'] are all equivalent.

A keyword (and its corresponding Card) can be deleted using the same index/name
syntax:

    >>> del prihdr[3] # delete the 2nd keyword
    >>> del prihdr['abc'] # get the value of the keyword 'abc'

Note that, like a regular Python list, the indexing updates after each delete,
so if ``del prihdr[3]`` is done two times in a row, the 2nd and 3rd keywords
are removed from the original header.

Slices are not accepted by the header attribute, so it is not possible to do
del ``prihdr[3:5]``, for example.

The method ``update(key, value, comment)`` is a more versatile way to update
keywords. It has the flexibility to update an existing keyword and in case the
keyword does not exist, add it to the header. It also allows the use to update
both the value and its comment. If it is a new keyword, the user can also
specify where to put it, using the before or after optional argument. The
default is to append at the end of the header.

    >>> prihdr.update('target', 'NGC1234', 'target name')
    >>> # place the next new keyword before the 'target' keyword
    >>> prihdr.update('newkey', 666, before='target') # comment is optional
    >>> # place the next new keyword after the 21st keyword
    >>> prihdr.update('newkey2', 42.0, 'another new key', after=20)


COMMENT, HISTORY, and Blank Keywords
------------------------------------

Most keywords in a FITS header have unique names. If there are more than two
cards sharing the same name, it is the first one accessed when referred by
name. The duplicates can only be accessed by numeric indexing.

There are three special keywords (their associated cards are sometimes referred
to as commentary cards), which commonly appear in FITS headers more than once.
They are (1) blank keyword, (2) HISTORY, and (3) COMMENT. Again, to get their
values (except for the first one), a user must use indexing.

The following header methods are provided in PyFITS to add new commentary
cards: `Header.add_history()`, `Header.add_comment()`, and
`Header.add_blank()`. They are provided because the `Header.update()` method
will not work - it will replace the first card of the same keyword.

Users can control where in the header to add the new commentary card(s) by
using the optional before and after arguments, similar to the ``update()``
method used for regular cards. If no before or after is specified, the new card
will be placed after the last one of the same kind (except  blank-key cards
which will always be placed at the end). If no card of the same kind exists, it
will be placed at the end. Here is an example:

    >>> hdu.header.add_history('history 1')
    >>> hdu.header.add_blank('blank 1')
    >>> hdu.header.add_comment('comment 1')
    >>> hdu.header.add_history('history 2')
    >>> hdu.header.add_blank('blank 2')
    >>> hdu.header.add_comment('comment 2'))

and the part in the modified header becomes:

.. parsed-literal::

    HISTORY history 1
    HISTORY history 2
            blank 1
    COMMENT comment 1
    COMMENT comment 2
            blank 2

Ironically, there is no comment in a commentary card , only a string value.


Card Images
===========

A FITS header consists of card images.

A card images in a FITS header consists of a keyword name, a value, and
optionally a comment. Physically, it takes 80 columns (bytes) - without
carriage return - in a FITS file's storage form. In PyFITS, each card image is
manifested by a Card object. There are also special kinds of cards: commentary
cards (see above) and card images taking more than one 80-column card image.
The latter will be discussed later.

Most of the time, a new Card object is created with the Card constructor:
``Card(key, value, comment)``. For example:

    >>> c1 = pyfits.Card('temp', 80.0, 'temperature, floating value')
    >>> c2 = pyfits.Card('detector', 1) # comment is optional
    >>> c3 = pyfits.Card('mir_revr', True, 'mirror reversed? Boolean value)
    >>> c4 = pyfits.Card('abc', 2+3j, 'complex value')
    >>> c5 = pyfits.Card('observer', 'Hubble', 'string value')

    >>> print c1; print c2; print c3; print c4; print c5 # show the card images
    TEMP = 80.0 / temperature, floating value
    DETECTOR= 1 /
    MIR_REVR= T / mirror reversed? Boolean value
    ABC = (2.0, 3.0) / complex value
    OBSERVER= 'Hubble ' / string value

Cards have the attributes ``.key``, ``.value``, and ``.comment``. Both
``.value`` and ``.comment`` can be changed but not the ``.key`` attribute.

The `Card()` constructor will check if the arguments given are conforming to
the FITS standard and has a fixed card image format. If the user wants to
create a card with a customized format or even a card which is not conforming
to the FITS standard (e.g. for testing purposes), the `Card.fromstring()`
method can be used.

Cards can be verified with `Card.verify()`. The non-standard card ``c2`` in the
example below, is flagged by such verification. More about verification in
PyFITS will be discussed in a later chapter.

    >>> c1 = pyfits.Card().fromstring('ABC = 3.456D023')
    >>> c2 = pyfits.Card().fromstring("P.I. ='Hubble'")
    >>> print c1; print c2
    ABC = 3.456D023
    P.I. ='Hubble'
    >>> c2.verify()
    Output verification result:
    Unfixable error: Illegal keyword name 'P.I.'


Card List
=========

The Header itself only has limited functionality. Many lower level operations
can only be achieved by going through its `CardList` object.

The header is basically a list of `Card` objects. This list can be manifested
as a `CardList` object in PyFITS. It is accessed via the `Header.ascardlist()`
method (or the ``.ascard`` attribute, for short). Since the header attribute
only refers to a card value, so when a user needs to access a card's other
properties (e.g. the comment) in a header, it has to go through the `CardList`.

Like the header's item, the `CardList`'s item can be accessed through either
the keyword name or index.

    >>> cards = prihdr.header.ascardlist()
    >>> cards['abc'].comment = 'new comment' # update the keyword ABC's comment
    >>> cards[3].key # see the keyword name of the 4th card
    >>> cards[10:20].keys() # see keyword names from cards 11 to 20


CONTINUE Cards
==============

The fact that the FITS standard only allows up to 8 characters for the keyword
name and 80 characters to contain the keyword, the value, and the comment is
restrictive for certain applications. To allow long string values for keywords,
a proposal was made in:

    http://legacy.gsfc.nasa.gov/docs/heasarc/ofwg/docs/ofwg_recomm/r13.html

by using the CONTINUE keyword after the regular 80-column containing the
keyword. PyFITS does support this convention, even though it is not a FITS
standard. The examples below show the use of CONTINUE is automatic for long
string values.

    >>> c = pyfits.Card('abc', 'abcdefg'*20)
    >>> print c
    ABC = 'abcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcd&'
    CONTINUE 'efgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefga&'
    CONTINUE 'bcdefg&'
    >>> c.value
    'abcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgab
    cdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefg'
    # both value and comments are long
    >>> c = pyfits.Card('abc', 'abcdefg'*10, 'abcdefg'*10)
    >>> print c
    ABC = 'abcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcd&'
    CONTINUE 'efg&'
    CONTINUE '&' / abcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefga
    CONTINUE '&' / bcdefg

Note that when CONTINUE card is used, at the end of each 80-characters card
image, an ampersand is present. The ampersand is not part of the string value.
Also, there is no "=" at the 9th column after CONTINUE. In the first example,
the entire 240 characters is considered a Card. So, if it is the nth card in a
header, the (n+1)th card refers to the next keyword, not the 80-characters
containing CONTINUE. These keywords having long string values can be accessed
and updated just like regular keywords.


HIERARCH Cards
==============

For keywords longer than 8 characters, there is a convention originated at ESO
to facilitate such use. It uses a special keyword HIERARCH with the actual long
keyword following. PyFITS supports this convention as well.

When creating or updating using the `Header.update()` method, it is necessary
to prepend 'hierarch' (case insensitive). But if the keyword is already in the
header, it can be accessed or updated by assignment by using the keyword name
diretly, with or without the 'hierarch' prepending.  The keyword name will
preserve its cases from its constructor, but when referring to the keyword, it
is case insensitive.

Examples follow:

    >>> c = pyfits.Card('abcdefghi', 10)
    ...
    ValueError: keyword name abcdefghi is too long (> 8), use HIERARCH.
    >>> c = pyfits.Card('hierarch abcdefghi', 10)
    >>> print c
    HIERARCH abcdefghi = 10
    >>> h = pyfits.PrimaryHDU()
    >>> h.header.update('hierarch abcdefghi', 99)
    >>> h.header.update('hierarch abcdefghi', 99)
    >>> h.header['abcdefghi']
    99
    >>> h.header['abcdefghi'] = 10
    >>> h.header['hierarch abcdefghi']
    10
    # case insensitive
    >>> h.header.update('hierarch ABCdefghi', 1000)
    >>> print h.header
    SIMPLE = T / conforms to FITS standard
    BITPIX = 8 / array data type
    NAXIS = 0 / number of array dimensions
    EXTEND = T
    HIERARCH ABCdefghi = 1000
    >>> h.header['hierarch abcdefghi']
    1000
