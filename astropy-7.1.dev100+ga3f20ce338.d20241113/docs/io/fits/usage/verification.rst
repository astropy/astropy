.. currentmodule:: astropy.io.fits

..  _fits_io_verification:

Verification
************

``astropy`` has built in a flexible scheme to verify FITS data conforming to
the FITS standard. The basic verification philosophy in ``astropy`` is to be
tolerant with input and strict with output.

When ``astropy`` reads a FITS file which does not conform to FITS standard, it
will not raise an error and exit. It will try to make the best educated
interpretation and only gives up when the offending data is accessed and no
unambiguous interpretation can be reached.

On the other hand, when writing to an output FITS file, the content to be
written must be strictly compliant to the FITS standard by default. This
default behavior can be overwritten by several other options, so the user will
not be held up because of a minor standard violation.


FITS Standard
=============

Since FITS standard is a "loose" standard, there are many places the violation
can occur and to enforce them all will be almost impossible. It is not uncommon
for major observatories to generate data products which are not 100% FITS
compliant. Some observatories have also developed their own nonstandard
dialect and some of these are so prevalent that they have become de facto
standards. Examples include the long string value and the use of the CONTINUE
card.

The violation of the standard can happen at different levels of the data
structure. ``astropy``'s verification scheme is developed on these hierarchical
levels. Here are the three ``astropy`` verification levels:

1. The HDU List

2. Each HDU

3. Each Card in the HDU Header

These three levels correspond to the three categories of objects:
:class:`HDUList`, any HDU (e.g., :class:`PrimaryHDU`, :class:`ImageHDU`, etc.),
and :class:`Card`. They are the only objects having the ``verify()`` method.
Most other classes in `astropy.io.fits` do not have a ``verify()`` method.

If ``verify()`` is called at the HDU List level, it verifies standard
compliance at all three levels, but a call of ``verify()`` at the Card level
will only check the compliance of that Card. Since ``astropy`` is tolerant when
reading a FITS file, no ``verify()`` is called on input. On output,
``verify()`` is called with the most restrictive option as the default.


Verification Options
====================

There are several options accepted by all verify(option) calls in ``astropy``.
In addition, they available for the ``output_verify`` argument of the following
methods: ``close()``, ``writeto()``, and ``flush()``. In these cases, they are
passed to a ``verify()`` call within these methods. The available options are:

**exception**

This option will raise an exception if any FITS standard is violated. This is
the default option for output (i.e., when ``writeto()``, ``close()``, or
``flush()`` is called). If a user wants to overwrite this default on output, the
other options listed below can be used.

**warn**

This option is the same as the ignore option but will send warning messages. It
will not try to fix any FITS standard violations whether fixable or not.

**ignore**

This option will ignore any FITS standard violation. On output, it will write
the HDU List content to the output FITS file, whether or not it is conforming
to the FITS standard.

The ignore option is useful in the following situations:

1. An input FITS file with nonstandard formatting is read and the user wants
   to copy or write out to an output file. The nonstandard formatting will be
   preserved in the output file.

2. A user wants to create a nonstandard FITS file on purpose, possibly for
   testing or consistency.

No warning message will be printed out. This is like a silent warning option
(see below).

**fix**

This option will try to fix any FITS standard violations. It is not always
possible to fix such violations. In general, there are two kinds of FITS
standard violations: fixable and non-fixable. For example, if a keyword has a
floating number with an exponential notation in lower case 'e' (e.g., 1.23e11)
instead of the upper case 'E' as required by the FITS standard, it is a fixable
violation. On the other hand, a keyword name like 'P.I.' is not fixable, since
it will not know what to use to replace the disallowed periods. If a violation
is fixable, this option will print out a message noting it is fixed. If it is
not fixable, it will throw an exception.

The principle behind fixing is to do no harm. For example, it is plausible to
'fix' a Card with a keyword name like 'P.I.' by deleting it, but ``astropy``
will not take such action to hurt the integrity of the data.

Not all fixes may be the "correct" fix, but at least ``astropy`` will try to
make the fix in such a way that it will not throw off other FITS readers.

**silentfix**

Same as fix, but will not print out informative messages. This may be useful in
a large script where the user does not want excessive harmless messages. If the
violation is not fixable, it will still throw an exception.

In addition the following combined options are available:

 * **fix+ignore**
 * **fix+warn**
 * **fix+exception**
 * **silentfix+ignore**
 * **silentfix+warn**
 * **silentfix+exception**

These options combine the semantics of the basic options. For example,
``silentfix+exception`` is actually equivalent to just ``silentfix`` in that
fixable errors will be fixed silently, but any unfixable errors will raise an
exception. On the other hand, ``silentfix+warn`` will issue warnings for
unfixable errors, but will stay silent about any fixed errors.


Verifications at Different Data Object Levels
=============================================

We will examine what ``astropy``'s verification does at the three different
levels:


Verification at HDUList
-----------------------

At the HDU List level, the verification is only for two simple cases:

1. Verify that the first HDU in the HDU list is a primary HDU. This is a
   fixable case. The fix is to insert a minimal primary HDU into the HDU list.

2. Verify the second or later HDU in the HDU list is not a primary HDU.
   Violation will not be fixable.


Verification at Each HDU
------------------------

For each HDU, the mandatory keywords, their locations in the header, and their
values will be verified. Each FITS HDU has a fixed set of required keywords in
a fixed order. For example, the primary HDU's header must at least have the
following keywords:

.. parsed-literal::

    SIMPLE =                     T /
    BITPIX =                     8 /
    NAXIS  =                     0

If any of the mandatory keywords are missing or in the wrong order, the fix
option will fix them::

    >>> from astropy.io import fits
    >>> filename = fits.util.get_testdata_filepath('verify.fits')
    >>> hdul = fits.open(filename)
    >>> hdul[0].header
    SIMPLE  =                    T / conforms to FITS standard
    NAXIS   =                    0 / NUMBER OF AXES
    BITPIX  =                    8 / BITS PER PIXEL
    >>> hdul[0].verify('fix') # doctest: +SHOW_WARNINGS
    VerifyWarning: Verification reported errors:
    VerifyWarning: 'BITPIX' card at the wrong place (card 2).
      Fixed by moving it to the right place (card 1).
    VerifyWarning: Note: astropy.io.fits uses zero-based indexing.
    >>> hdul[0].header           # voila!
    SIMPLE  =                    T / conforms to FITS standard
    BITPIX  =                    8 / BITS PER PIXEL
    NAXIS   =                    0 / NUMBER OF AXES
    >>> hdul.close()

Verification at Each Card
-------------------------

The lowest level, the Card, also has the most complicated verification
possibilities.

Examples
^^^^^^^^

..
  EXAMPLE START
  Verification at Each Card in astropy.io.fits

Here is a list of fixable and not fixable Cards:

Fixable Cards:

1. Floating point numbers with lower case 'e' or 'd'::

    >>> from astropy.io import fits
    >>> c = fits.Card.fromstring('FIX1    = 2.1e23')
    >>> c.verify('silentfix')
    >>> print(c)
    FIX1    =               2.1E23

2. The equal sign is before column nine in the card image::

    >>> c = fits.Card.fromstring('FIX2= 2')
    >>> c.verify('silentfix')
    >>> print(c)
    FIX2    =                    2

3. String value without enclosing quotes::

    >>> c = fits.Card.fromstring('FIX3    = string value without quotes')
    >>> c.verify('silentfix')
    >>> print(c)
    FIX3    = 'string value without quotes'

4. Missing equal sign before column nine in the card image.

5. Space between numbers and E or D in floating point values::

    >>> c = fits.Card.fromstring('FIX5    = 2.4 e 03')
    >>> c.verify('silentfix')
    >>> print(c)
    FIX5    =               2.4E03

6. Unparsable values will be "fixed" as a string::

    >>> c = fits.Card.fromstring('FIX6    = 2 10 ')
    >>> c.verify('fix+warn') # doctest: +SHOW_WARNINGS
    VerifyWarning: Verification reported errors:
    VerifyWarning: Card 'FIX6' is not FITS standard
     (invalid value string: '2 10').
       Fixed 'FIX6' card to meet the FITS standard.
    VerifyWarning: Note: astropy.io.fits uses zero-based indexing.
    >>> print(c)
    FIX6    = '2 10    '

Unfixable Cards:

1. Illegal characters in keyword name.

We will summarize the verification with a "life-cycle" example::

    >>> h = fits.PrimaryHDU()  # create a PrimaryHDU
    >>> # Try to add an non-standard FITS keyword 'P.I.' (FITS does no allow
    >>> # '.' in the keyword), if using the update() method - doesn't work!
    >>> h.header['P.I.'] = 'Hubble' # doctest: +SHOW_WARNINGS
    VerifyWarning: Keyword name 'P.I.' is greater than 8 characters or
     contains characters not allowed by the FITS standard;
      a HIERARCH card will be created.
    >>> # Have to do it the hard way (so a user will not do this by accident)
    >>> # First, create a card image and give verbatim card content (including
    >>> # the proper spacing, but no need to add the trailing blanks)
    >>> c = fits.Card.fromstring("P.I. = 'Hubble'")
    >>> h.header.append(c)  # then append it to the header
    >>> # Now if we try to write to a FITS file, the default output
    >>> # verification will not take it.
    >>> h.writeto('pi.fits')  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
     ...
    VerifyError: HDU 0:
        Card 5:
            Card 'P.I. ' is not FITS standard (equal sign not at column 8).
            Illegal keyword name 'P.I. '
    >>> # Must set the output_verify argument to 'ignore', to force writing a
    >>> # non-standard FITS file
    >>> h.writeto('pi.fits', output_verify='ignore')
    >>> # Now reading a non-standard FITS file
    >>> # astropy.io.fits is magnanimous in reading non-standard FITS files
    >>> hdul = fits.open('pi.fits')
    >>> hdul[0].header # doctest: +SHOW_WARNINGS
    SIMPLE  =            T / conforms to FITS standard
    BITPIX  =            8 / array data type
    NAXIS   =            0 / number of array dimensions
    EXTEND  =            T
    HIERARCH P.I. = 'Hubble  '
    P.I.    = 'Hubble  '
    VerifyWarning: Verification reported errors:
    VerifyWarning: Card 'P.I. ' is not FITS standard (equal sign
     not at column 8).  Fixed 'P.I. ' card to meet the FITS standard.
    VerifyWarning: Unfixable error: Illegal keyword name 'P.I. '
    VerifyWarning: Note: astropy.io.fits uses zero-based indexing.
    >>> # even when you try to access the offending keyword, it does NOT
    >>> # complain
    >>> hdul[0].header['p.i.']
    'Hubble'
    >>> # But if you want to make sure if there is anything wrong/non-standard,
    >>> # use the verify() method
    >>> hdul.verify() # doctest: +SHOW_WARNINGS
    VerifyWarning: Verification reported errors:
    VerifyWarning: HDU 0:
    VerifyWarning:     Card 5:
    VerifyWarning:         Illegal keyword name 'P.I. '
    VerifyWarning: Note: astropy.io.fits uses zero-based indexing.
    >>> hdul.close()

..
  EXAMPLE END

Verification Using the FITS Checksum Keyword Convention
=======================================================

The North American FITS committee has reviewed the FITS Checksum Keyword
Convention for possible adoption as a FITS Standard. This convention provides
an integrity check on information contained in FITS HDUs. The convention
consists of two header keyword cards: CHECKSUM and DATASUM. The CHECKSUM
keyword is defined as an ASCII character string whose value forces the 32-bit
1's complement checksum accumulated over all the 2880-byte FITS logical records
in the HDU to equal negative zero. The DATASUM keyword is defined as a
character string containing the unsigned integer value of the 32-bit 1's
complement checksum of the data records in the HDU. Verifying the
accumulated checksum is still equal to negative zero provides a fairly reliable
way to determine that the HDU has not been modified by subsequent data
processing operations or corrupted while copying or storing the file on
physical media.

In order to avoid any impact on performance, by default ``astropy`` will not
verify HDU checksums when a file is opened or generate checksum values when a
file is written. In fact, CHECKSUM and DATASUM cards are automatically removed
from HDU headers when a file is opened, and any CHECKSUM or DATASUM cards are
stripped from headers when an HDU is written to a file. In order to verify the
checksum values for HDUs when opening a file, the user must supply the checksum
keyword argument in the call to the open convenience function with a value of
True. When this is done, any checksum verification failure will cause a
warning to be issued (via the warnings module). If checksum verification is
requested in the open, and no CHECKSUM or DATASUM cards exist in the HDU
header, the file will open without comment. Similarly, in order to output the
CHECKSUM and DATASUM cards in an HDU header when writing to a file, the user
must supply the checksum keyword argument with a value of True in the call to
the ``writeto()`` function. It is possible to write only the DATASUM card to the
header by supplying the checksum keyword argument with a value of 'datasum'.

Examples
--------

..
  EXAMPLE START
  Verification Using the FITS Checksum Keyword Convention

To verify the checksum values for HDUs when opening a file::

    >>> # Open the file checksum.fits verifying the checksum values for all HDUs
    >>> filename = fits.util.get_testdata_filepath('checksum.fits')
    >>> hdul = fits.open(filename, checksum=True)
    >>> hdul.close()
    >>> # Open the file in.fits where checksum verification fails
    >>> filename = fits.util.get_testdata_filepath('checksum_false.fits')
    >>> hdul = fits.open(filename, checksum=True) # doctest: +SHOW_WARNINGS
    AstropyUserWarning: Checksum verification failed for HDU ('PRIMARY', 1).
    AstropyUserWarning: Datasum verification failed for HDU ('PRIMARY', 1).
    AstropyUserWarning: Checksum verification failed for HDU ('RATE', 1).
    AstropyUserWarning: Datasum verification failed for HDU ('RATE', 1).
    >>> # Create file out.fits containing an HDU constructed from data
    >>> # containing both CHECKSUM and DATASUM cards.
    >>> data = hdul[0].data
    >>> fits.writeto('out.fits', data=data, checksum=True)
    >>> hdun = fits.open('out.fits', checksum=True)
    >>> hdun.close()

    >>> # Create file out.fits containing all the HDUs in the HDULIST
    >>> # hdul with each HDU header containing only the DATASUM card
    >>> hdul.writeto('out2.fits', checksum='datasum')

    >>> # Create file out.fits containing the HDU hdu with both CHECKSUM
    >>> # and DATASUM cards in the header
    >>> hdu = hdul[1]
    >>> hdu.writeto('out3.fits', checksum=True)

    >>> # Append a new HDU constructed from array data to the end of
    >>> # the file existingfile.fits with only the appended HDU
    >>> # containing both CHECKSUM and DATASUM cards.
    >>> fits.append('out3.fits', data, checksum=True)
    >>> hdul.close()

..
  EXAMPLE END
