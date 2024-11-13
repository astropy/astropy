.. currentmodule:: astropy.io.fits

FITS Headers
************

In the next three chapters, more detailed information including examples will
be explained for manipulating FITS headers, image/array data, and table data
respectively.


Header of an HDU
================

Every Header Data Unit (HDU) normally has two components: header and data. In
``astropy`` these two components are accessed through the two attributes of the
HDU, ``hdu.header`` and ``hdu.data``.

While an HDU may have empty data (i.e., the ``.data`` attribute is `None`), any
HDU will always have a header. When an HDU is created with a constructor (e.g.,
``hdu = PrimaryHDU(data, header)``), the user may supply the header value from
an existing HDU's header and the data value from a ``numpy`` array. If the
defaults (None) are used, the new HDU will have the minimal required keywords
for an HDU of that type::

    >>> from astropy.io import fits
    >>> hdu = fits.PrimaryHDU()
    >>> hdu.header  # show the all of the header cards
    SIMPLE  =                    T / conforms to FITS standard
    BITPIX  =                    8 / array data type
    NAXIS   =                    0 / number of array dimensions
    EXTEND  =                    T

A user can use any header and any data to construct a new HDU. ``astropy`` will
strip any keywords that describe the data structure leaving only your
informational keywords. Later it will add back in the required structural
keywords for compatibility with the new HDU and any data added to it. So, a
user can use a table HDU's header to construct an image HDU and vice versa. The
constructor will also ensure the data type and dimension information in the
header agree with the data.


The Header Attribute
====================

Value Access, Updating, and Creating
------------------------------------

As shown in the :ref:`Getting Started <tutorial>` tutorial, keyword values can
be accessed via keyword name or index of an HDU's header attribute. You can
also use the wildcard character ``*`` to get the keyword value pairs that match
your search string. Here is a quick summary::

    >>> fits_image_filename = fits.util.get_testdata_filepath('test0.fits')
    >>> hdul = fits.open(fits_image_filename)  # open a FITS file
    >>> hdr = hdul[0].header  # the primary HDU header
    >>> print(hdr[34])  # get the 2nd keyword's value
    96
    >>> hdr[34] = 20  # change its value
    >>> hdr['DARKCORR']  # get the value of the keyword 'darkcorr'
    'OMIT'
    >>> hdr['DARKCOR*']  # get keyword values using wildcard matching
    DARKCORR= 'OMIT              ' / Do dark correction: PERFORM, OMIT, COMPLETE
    >>> hdr['darkcorr'] = 'PERFORM'  # change darkcorr's value

Keyword names are case-insensitive except in a few special cases (see the
sections on HIERARCH card and record-valued cards). Thus, ``hdr['abc']``,
``hdr['ABC']``, or ``hdr['aBc']`` are all equivalent.

As with Python's :class:`dict` type, new keywords can also be added to the
header using assignment syntax::

    >>> hdr = hdul[1].header
    >>> 'DARKCORR' in hdr  # Check for existence
    False
    >>> hdr['DARKCORR'] = 'OMIT'  # Add a new DARKCORR keyword

You can also add a new value *and* comment by assigning them as a tuple::

    >>> hdr['DARKCORR'] = ('OMIT', 'Dark Image Subtraction')

.. note::

    An important point to note when adding new keywords to a header is that by
    default they are not appended *immediately* to the end of the file.
    Rather, they are appended to the last non-commentary keyword. This is in
    order to support the common use case of always having all HISTORY keywords
    grouped together at the end of a header. A new non-commentary keyword will
    be added at the end of the existing keywords, but before any
    HISTORY/COMMENT keywords at the end of the header.

    There are a couple of ways to override this functionality:

    * Use the :meth:`Header.append` method with the ``end=True`` argument:

        >>> hdr.append(('DARKCORR', 'OMIT', 'Dark Image Subtraction'), end=True)

      This forces the new keyword to be added at the actual end of the header.

    * The :meth:`Header.insert` method will always insert a new keyword exactly
      where you ask for it:

        >>> del hdr['DARKCORR']  # Delete previous insertion for doctest
        >>> hdr.insert(20, ('DARKCORR', 'OMIT', 'Dark Image Subtraction'))

      This inserts the DARKCORR keyword before the 20th keyword in the
      header no matter what it is.

A keyword (and its corresponding card) can be deleted using the same index/name
syntax::

    >>> del hdr[3]  # delete the 2nd keyword
    >>> del hdr['DARKCORR']  # delete the value of the keyword 'DARKCORR'

Note that, like a regular Python list, the indexing updates after each delete,
so if ``del hdr[3]`` is done two times in a row, the fourth and fifth keywords
are removed from the original header. Likewise, ``del hdr[-1]`` will delete
the last card in the header.

It is also possible to delete an entire range of cards using the slice syntax::

    >>> del hdr[3:5]

The method :meth:`Header.set` is another way to update the value or comment
associated with an existing keyword, or to create a new keyword. Most of its
functionality can be duplicated with the dict-like syntax shown above. But in
some cases it might be more clear. It also has the advantage of allowing a user
to either move cards within the header or specify the location of a new card
relative to existing cards::

    >>> hdr.set('target', 'NGC1234', 'target name')
    >>> # place the next new keyword before the 'TARGET' keyword
    >>> hdr.set('newkey', 666, before='TARGET')  # comment is optional
    >>> # place the next new keyword after the 21st keyword
    >>> hdr.set('newkey2', 42.0, 'another new key', after=20)

In FITS headers, each keyword may also have a comment associated with it
explaining its purpose. The comments associated with each keyword are accessed
through the :attr:`~Header.comments` attribute::

    >>> hdr['NAXIS']
    2
    >>> hdr.comments['NAXIS']
    'number of data axes'
    >>> hdr.comments['NAXIS'] = 'The number of image axes'  # Update
    >>> hdul.close()  # close the HDUList again

Comments can be accessed in all of the same ways that values are accessed,
whether by keyword name or card index. Slices are also possible. The only
difference is that you go through ``hdr.comments`` instead of just ``hdr`` by
itself.


COMMENT, HISTORY, and Blank Keywords
------------------------------------

Most keywords in a FITS header have unique names. If there are more than two
cards sharing the same name, it is the first one accessed when referred by
name. The duplicates can only be accessed by numeric indexing.

There are three special keywords (their associated cards are sometimes referred
to as commentary cards), which commonly appear in FITS headers more than once.
They are (1) blank keyword, (2) HISTORY, and (3) COMMENT. Unlike other
keywords, when accessing these keywords they are returned as a list::

    >>> filename = fits.util.get_testdata_filepath('history_header.fits')
    >>> with fits.open(filename) as hdul:  # open a FITS file
    ...     hdr = hdul[0].header

    >>> hdr['HISTORY']
    I updated this file on 02/03/2011
    I updated this file on 02/04/2011

These lists can be sliced like any other list. For example, to display just the
last HISTORY entry, use ``hdr['history'][-1]``. Existing commentary cards
can also be updated by using the appropriate index number for that card.

New commentary cards can be added like any other card by using the dict-like
keyword assignment syntax, or by using the :meth:`Header.set` method. However,
unlike with other keywords, a new commentary card is always added and appended
to the last commentary card with the same keyword, rather than to the end of
the header.

Example
^^^^^^^

..
  EXAMPLE START
  Manipulating FITS Headers in astropy.io.fits

To add a new commentary card::

    >>> hdu.header['HISTORY'] = 'history 1'
    >>> hdu.header[''] = 'blank 1'
    >>> hdu.header['COMMENT'] = 'comment 1'
    >>> hdu.header['HISTORY'] = 'history 2'
    >>> hdu.header[''] = 'blank 2'
    >>> hdu.header['COMMENT'] = 'comment 2'

and the part in the modified header becomes:

.. parsed-literal::

    HISTORY history 1
    HISTORY history 2
            blank 1
            blank 2
    COMMENT comment 1
    COMMENT comment 2


Users can also directly control exactly where in the header to add a new
commentary card by using the :meth:`Header.insert` method.

.. note::

    Ironically, there is no comment in a commentary card, only a string
    value.

..
  EXAMPLE END

Undefined Values
----------------

FITS headers can have undefined values and these are represented in Python
with the special value `None`. `None` can be used when assigning values
to a `~astropy.io.fits.Header` or `~astropy.io.fits.Card`.

    >>> hdr = fits.Header()
    >>> hdr['UNDEF'] = None
    >>> hdr['UNDEF'] is None
    True
    >>> repr(hdr)
    'UNDEF   =                                                                       '
    >>> hdr.append('UNDEF2')
    >>> hdr['UNDEF2'] is None
    True
    >>> hdr.append(('UNDEF3', None, 'Undefined value'))
    >>> str(hdr.cards[-1])
    'UNDEF3  =  / Undefined value                                                    '


Card Images
===========

A FITS header consists of card images.

A card image in a FITS header consists of a keyword name, a value, and
optionally a comment. Physically, it takes 80 columns (bytes) — without carriage
return — in a FITS file's storage format. In ``astropy``, each card image is
manifested by a :class:`Card` object. There are also special kinds of cards:
commentary cards (see above) and card images taking more than one 80-column
card image. The latter will be discussed later.

Most of the time the details of dealing with cards are handled by the
:class:`Header` object, and it is not necessary to directly manipulate cards.
In fact, most :class:`Header` methods that accept a ``(keyword, value)`` or
``(keyword, value, comment)`` tuple as an argument can also take a
:class:`Card` object as an argument. :class:`Card` objects are just wrappers
around such tuples that provide the logic for parsing and formatting individual
cards in a header. There is usually nothing gained by manually using a
:class:`Card` object, except to examine how a card might appear in a header
before actually adding it to the header.

A new Card object is created with the :class:`Card` constructor:
``Card(key, value, comment)``.

Example
-------

..
  EXAMPLE START
  Card Images in FITS Headers in astropy.io.fits

To create a new Card object::

    >>> c1 = fits.Card('TEMP', 80.0, 'temperature, floating value')
    >>> c2 = fits.Card('DETECTOR', 1)  # comment is optional
    >>> c3 = fits.Card('MIR_REVR', True,
    ...                'mirror reversed? Boolean value')
    >>> c4 = fits.Card('ABC', 2+3j, 'complex value')
    >>> c5 = fits.Card('OBSERVER', 'Hubble', 'string value')

    >>> print(c1); print(c2); print(c3); print(c4); print(c5)  # show the cards
    TEMP    =                 80.0 / temperature, floating value
    DETECTOR=                    1
    MIR_REVR=                    T / mirror reversed? Boolean value
    ABC     =           (2.0, 3.0) / complex value
    OBSERVER= 'Hubble  '           / string value

Cards have the attributes ``.keyword``, ``.value``, and ``.comment``. Both
``.value`` and ``.comment`` can be changed but not the ``.keyword`` attribute.
In other words, once a card is created, it is created for a specific, immutable
keyword.

The :meth:`Card` constructor will check if the arguments given are conforming
to the FITS standard and has a fixed card image format. If the user wants to
create a card with a customized format or even a card which is not conforming
to the FITS standard (e.g., for testing purposes), the :meth:`Card.fromstring`
class method can be used.

Cards can be verified with :meth:`Card.verify`. The nonstandard card ``c2`` in
the example below is flagged by such verification. More about verification in
``astropy`` will be discussed in a later chapter.

::

    >>> c1 = fits.Card.fromstring('ABC     = 3.456D023')
    >>> c2 = fits.Card.fromstring("P.I. ='Hubble'")
    >>> print(c1)
    ABC     = 3.456D023
    >>> print(c2)  # doctest: +SKIP
    P.I. ='Hubble'
    >>> c2.verify()  # doctest: +SKIP
    Output verification result:
    Unfixable error: Illegal keyword name 'P.I.'

A list of the :class:`Card` objects underlying a :class:`Header` object can be
accessed with the :attr:`Header.cards` attribute. This list is only meant for
observing, and should not be directly manipulated. In fact, it is only a
copy — modifications to it will not affect the header from which it came. Use
the methods provided by the :class:`Header` class instead.

..
  EXAMPLE END


CONTINUE Cards
==============

The fact that the FITS standard only allows up to eight characters for the
keyword name and 80 characters to contain the keyword, the value, and the
comment is restrictive for certain applications. To allow long string values
for keywords, a proposal was made in:

    https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/ofwg_recomm/r13.html

by using the CONTINUE keyword after the regular 80 column containing the
keyword. ``astropy`` does support this convention, which is a part of the FITS
standard since version 4.0.

Examples
--------

..
  EXAMPLE START
  CONTINUE Cards for Long String Values in astropy.io.fits

The examples below show that the use of CONTINUE is automatic for long
string values::

    >>> hdr = fits.Header()
    >>> hdr['abc'] = 'abcdefg' * 20
    >>> hdr
    ABC     = 'abcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcd&'
    CONTINUE  'efgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefga&'
    CONTINUE  'bcdefg'
    >>> hdr['abc']
    'abcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefg'
    >>> # both value and comments are long
    >>> hdr['abc'] = ('abcdefg' * 10, 'abcdefg' * 10)
    >>> hdr
    ABC     = 'abcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcd&'
    CONTINUE  'efg&'
    CONTINUE  '&' / abcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefgabcdefga
    CONTINUE  '' / bcdefg

Note that when a CONTINUE card is used, at the end of each 80-character card
image, an ampersand is present. The ampersand is not part of the string value.
Also, there is no "=" at the ninth column after CONTINUE. In the first example,
the entire 240 characters is treated by ``astropy`` as a single card. So, if it
is the nth card in a header, the (n+1)th card refers to the next keyword, not
the next CONTINUE card. As such, CONTINUE cards are transparently handled by
``astropy`` as a single logical card, and it is generally not necessary to worry
about the details of the format. Keywords that resolve to a set of CONTINUE
cards can be accessed and updated just like regular keywords.

..
  EXAMPLE END


HIERARCH Cards
==============

For keywords longer than eight characters, there is a convention originated at
the European Southern Observatory (ESO) to facilitate such use. It uses a
special keyword HIERARCH with the actual long keyword following. ``astropy``
supports this convention as well.

If a keyword contains more than eight characters ``astropy`` will automatically
use a HIERARCH card, but will also issue a warning in case this is in error.
However, you may explicitly request a HIERARCH card by prepending the keyword
with 'HIERARCH ' (just as it would appear in the header). For example,
``hdr['HIERARCH abcdefghi']`` will create the keyword ``abcdefghi`` without
displaying a warning. Once created, HIERARCH keywords can be accessed like any
other: ``hdr['abcdefghi']``, without prepending 'HIERARCH' to the keyword.

Examples
--------

..
  EXAMPLE START
  HIERARCH Cards for Keywords Longer than Eight Characters in astropy.io.fits

``astropy`` will use a HIERARCH card and issue a warning when keywords contain
more than eight characters::

    >>> # this will result in a Warning because a HIERARCH card is implicitly created
    >>> c = fits.Card('abcdefghi', 10)  # doctest: +SKIP
    >>> print(c)  # doctest: +SKIP
    HIERARCH abcdefghi = 10
    >>> c = fits.Card('hierarch abcdefghi', 10)
    >>> print(c)
    HIERARCH abcdefghi = 10
    >>> hdu = fits.PrimaryHDU()
    >>> hdu.header['hierarch abcdefghi'] =  99
    >>> hdu.header['abcdefghi']
    99
    >>> hdu.header['abcdefghi'] = 10
    >>> hdu.header['abcdefghi']
    10
    >>> hdu.header
    SIMPLE  =                    T / conforms to FITS standard
    BITPIX  =                    8 / array data type
    NAXIS   =                    0 / number of array dimensions
    EXTEND  =                    T
    HIERARCH abcdefghi = 10

..
  EXAMPLE END

.. note::

    A final point to keep in mind about the :class:`Header` class is that much
    of its design is intended to abstract away quirks about the FITS format.
    This is why, for example, it will automatically create CONTINUE and
    HIERARCH cards. The Header is just a data structure, and as a user you
    should not have to worry about how it ultimately gets serialized to a header
    in a FITS file.

    Though there are some areas where it is almost impossible to hide away the
    quirks of the FITS format, ``astropy`` tries to make it so that you have to
    think about it as little as possible. If there are any areas that are left
    vague or difficult to understand about how the header is constructed, please
    let us know, as there are probably areas where this can be
    improved on even more.
