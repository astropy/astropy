.. currentmodule:: astropy.io.fits
.. doctest-skip-all

.. _header-transition-guide:

*********************************
Header Interface Transition Guide
*********************************

.. note::

    This guide was originally included with the release of PyFITS 3.1, and
    still references PyFITS in many places, though the examples have been
    updated for ``astropy.io.fits``. It is still useful here for informational
    purposes, though Astropy has always used the PyFITS 3.1 Header interface.

PyFITS v3.1 included an almost complete rewrite of the :class:`Header`
interface. Although the new interface is largely compatible with the old
interface (whether due to similarities in the design, or backwards-compatibility
support), there are enough differences that a full explanation of the new
interface is merited.

Background
==========

Prior to 3.1, PyFITS users interacted with FITS headers by way of three
different classes: :class:`Card`, ``CardList``, and :class:`Header`.

The Card class represents a single header card with a keyword, value, and
comment. It also contains all of the machinery for parsing FITS header cards,
given the 80-character string, or "card image" read from the header.

The CardList class is actually a subclass of Python's `list` built-in. It was
meant to represent the actual list of cards that make up a header. That is, it
represents an ordered list of cards in the physical order that they appear in
the header. It supports the usual list methods for inserting and appending new
cards into the list. It also supports `dict`-like keyword access, where
``cardlist['KEYWORD']`` would return the first card in the list with the given
keyword.

A lot of the functionality for manipulating headers was actually buried in the
CardList class. The Header class was more of a wrapper around CardList that
added a little bit of abstraction. It also implemented a partial dict-like
interface, though for Headers a keyword lookup returned the header value
associated with that keyword, not the Card object, and almost every
method on the Header class was just performing some operations on the
underlying CardList.

The problem was that there were certain things a user could *only* do by
directly accessing the CardList, such as look up the comments on a card or
access cards that have duplicate keywords, such as HISTORY. Another long-
standing misfeature was that slicing a Header object actually returned a
CardList object, rather than a new Header. For all but the simplest use cases,
working with CardList objects was largely unavoidable.

But it was realized that CardList is really an implementation detail
not representing any element of a FITS file distinct from the header itself.
Users familiar with the FITS format know what a header is, but it is not clear
how a "card list" is distinct from that, or why operations go through the
Header object, while some have to be performed through the CardList.

So the primary goal of this redesign was to eliminate the ``CardList`` class
altogether, and make it possible for users to perform all header manipulations
directly through :class:`Header` objects. It also tried to present headers as
similarly as possible to a more familiar data structure — an ordered mapping
(or :class:`~collections.OrderedDict` in Python) for ease of use by new users
less familiar with the FITS format, though there are still many added
complexities for dealing with the idiosyncrasies of the FITS format.


Deprecation Warnings
====================

A few older methods on the :class:`Header` class have been marked as deprecated,
either because they have been renamed to a more `PEP 8`_-compliant name, or
because have become redundant due to new features. To check if your code is
using any deprecated methods or features, run your code with ``python -Wd``.
This will output any deprecation warnings to the console.

Two of the most common deprecation warnings related to Headers are:

- ``Header.has_key``: this has been deprecated since PyFITS 3.0,
  just as Python's `dict.has_key` is deprecated. To check a key's presence
  in a mapping object like `dict` or :class:`Header`, use the ``key in d``
  syntax. This has long been the preference in Python.

- ``Header.ascardlist`` and ``Header.ascard``: these were used to
  access the ``CardList`` object underlying a header. They should still
  work, and return a skeleton CardList implementation that should support most
  of the old CardList functionality. But try removing as much of this as
  possible. If direct access to the :class:`Card` objects making up a header
  is necessary, use :attr:`Header.cards`, which returns an iterator over the
  cards. More on that below.

.. _PEP 8: https://www.python.org/dev/peps/pep-0008/

New Header Design
=================

The new :class:`Header` class is designed to work as a drop-in replacement for
a `dict` via `duck typing`_. That is, although it is not a subclass of `dict`,
it implements all of the same methods and interfaces. In particular, it is
similar to an :class:`~collections.OrderedDict` in that the order of insertions
is preserved. However, Header also supports many additional features and
behaviors specific to the FITS format. It should also be noted that while the
old Header implementation also had a dict-like interface, it did not implement
the *entire* dict interface as the new Header does.

Although the new Header is used like a dict/mapping in most cases, it also
supports a `list` interface. The list-like interface is a bit idiosyncratic in
that in some contexts the Header acts like a list of values, in others like a
list of keywords, and in a few contexts like a list of :class:`Card` objects.
This may be the most difficult aspect of the new design, but there is a logic
to it.

As with the old Header implementation, integer index access is supported:
``header[0]`` returns the value of the first keyword. However, the
:meth:`Header.index` method treats the header as though it is a list of
keywords and returns the index of a given keyword. For example::

    >>> header.index('BITPIX')
    2

:meth:`Header.count` is similar to `list.count` and also takes a keyword as
its argument::

    >>> header.count('HISTORY')
    20

A good rule of thumb is that any item access using square brackets ``[]``
returns *value* in the header, whether using keyword or index lookup. Methods
like :meth:`~Header.index` and :meth:`~Header.count` that deal with the order
and quantity of items in the Header generally work on keywords. Finally,
methods like :meth:`~Header.insert` and :meth:`~Header.append` that add new
items to the header work on cards.

Aside from the list-like methods, the new Header class works very similarly to
the old implementation for most basic use cases and should not present too many
surprises. There are differences, however:

- As before, the Header() initializer can take a list of :class:`Card` objects
  with which to fill the header. However, now any iterable may be used. It is
  also important to note that *any* Header method that accepts :class:`Card`
  objects can also accept 2-tuples or 3-tuples in place of Cards. That is,
  either a ``(keyword, value, comment)`` tuple or a ``(keyword, value)`` tuple
  (comment is assumed blank) may be used anywhere in place of a Card object.
  This is even preferred, as it involves less typing. For example::

      >>> from astropy.io import fits
      >>> header = fits.Header([('A', 1), ('B', 2), ('C', 3, 'A comment')])
      >>> header
      A       =                    1
      B       =                    2
      C       =                    3 / A comment

- As demonstrated in the previous example, the ``repr()`` for a Header (that is,
  the text that is displayed when entering a Header object in the Python
  console as an expression), shows the header as it would appear in a FITS file.
  This inserts newlines after each card so that it is readable regardless of
  terminal width. It is *not* necessary to use ``print header`` to view this.
  Entering ``header`` displays the header contents as it would appear in the
  file (sans the END card).

- ``len(header)`` is now supported (previously it was necessary to do
  ``len(header.ascard)``). This returns the total number of cards in the
  header, including blank cards, but excluding the END card.

- FITS supports having duplicate keywords, although they are generally in error
  except for commentary keywords like COMMENT and HISTORY. PyFITS now supports
  reading, updating, and deleting duplicate keywords; instead of using the
  keyword by itself, use a ``(keyword, index)`` tuple. For example,
  ``('HISTORY', 0)`` represents the first HISTORY card, ``('HISTORY', 1)``
  represents the second HISTORY card, and so on. In fact, when a keyword is
  used by itself, it is shorthand for ``(keyword, 0)``. It is now possible to
  delete an accidental duplicate like so::

      >>> del header[('NAXIS', 1)]

  This will remove an accidental duplicate NAXIS card from the header.

- Even if there are duplicate keywords, keyword lookups like
  ``header['NAXIS']`` will always return the value associated with the first
  copy of that keyword, with one exception: commentary keywords like COMMENT
  and HISTORY are expected to have duplicates. So ``header['HISTORY']``, for
  example, returns the whole sequence of HISTORY values in the correct order.
  This list of values can be sliced arbitrarily. For example, to view the last
  three history entries in a header::

      >>> hdulist[0].header['HISTORY'][-3:]
        reference table oref$laf13367o_pct.fits
        reference table oref$laf13369o_apt.fits
      Heliocentric correction = 16.225 km/s

- Subscript assignment can now be used to add new keywords to the header. Just
  as with a normal `dict`, ``header['NAXIS'] = 1`` will either update the NAXIS
  keyword if it already exists, or add a new NAXIS keyword with a value of
  ``1`` if it does not exist. In the old interface this would return a
  `KeyError` if NAXIS did not exist, and the only way to add a new
  keyword was through the update() method.

  By default, new keywords added in this manner are added to the end of the
  header, with a few FITS-specific exceptions:

  * If the header contains extra blank cards at the end, new keywords are added
    before the blanks.

  * If the header ends with a list of commentary cards — for example, a sequence
    of HISTORY cards — those are kept at the end, and new keywords are inserted
    before the commentary cards.

  * If the keyword is a commentary keyword like COMMENT or HISTORY (or an empty
    string for blank keywords), a *new* commentary keyword is always added and
    appended to the last commentary keyword of the same type. For example,
    HISTORY keywords are always placed after the last history keyword::

        >>> header = fits.Header()
        >>> header['COMMENT'] = 'Comment 1'
        >>> header['HISTORY'] = 'History 1'
        >>> header['COMMENT'] = 'Comment 2'
        >>> header['HISTORY'] = 'History 2'
        >>> header
        COMMENT Comment 1
        COMMENT Comment 2
        HISTORY History 1
        HISTORY History 2

  These behaviors represent a sensible default behavior for keyword assignment,
  and the same behavior as :meth:`~Header.update` in the old Header
  implementation. The default behaviors may still be bypassed through the use
  of other assignment methods like the :meth:`Header.set` and
  :meth:`Header.append` methods described later.

- It is now also possible to assign a value and a comment to a keyword
  simultaneously using a tuple::

      >>> header['NAXIS'] = (2, 'Number of axis')

  This will update the value and comment of an existing keyword, or add a new
  keyword with the given value and comment.

- There is a new :attr:`Header.comments` attribute which lists all of the
  comments associated with keywords in the header (not to be confused with
  COMMENT cards). This allows viewing and updating the comments on specific
  cards::

      >>> header.comments['NAXIS']
      Number of axis
      >>> header.comments['NAXIS'] = 'Number of axes'
      >>> header.comments['NAXIS']
      Number of axes

- When deleting a keyword from a header, do not assume that the keyword already
  exists. In the old Header implementation, this action would silently do
  nothing. For backwards-compatibility, it is still okay to delete a
  nonexistent keyword, but a warning will be raised. In the future this
  *will* be changed so that trying to delete a nonexistent keyword raises a
  `KeyError`. This is for consistency with the behavior of Python dicts. So
  unless you know for certain that a keyword exists before deleting it, it is
  best to do something like::

      >>> try:
      ...     del header['BITPIX']
      ... except KeyError:
      ...     pass

  Or if you prefer to look before you leap::

      >>> if 'BITPIX' in header:
      ...     del header['BITPIX']

- ``del header`` now supports slices. For example, to delete the last three
  keywords from a header::

      >>> del header[-3:]

- Two headers can now be compared for equality — previously no two Header
  objects were the same. Now they compare as equal if they contain the exact
  same content. That is, this requires strict equality.

- Two headers can now be added with the '+' operator, which returns a copy of
  the left header extended by the right header with :meth:`~Header.extend`.
  Assignment addition is also possible.

- The Header.update() method used commonly with the old Header API has been
  renamed to :meth:`Header.set`. The primary reason for this change is very
  simple: Header implements the `dict` interface, which already has a method
  called update(), but that behaves differently from the old Header.update().

  The details of the new update() can be read in the API docs, but it is very
  similar to `dict.update`. It also supports backwards compatibility with the
  old update() by analysis of the arguments passed to it, so existing code will
  not break immediately. However, this *will* cause a deprecation warning to
  be output if they are enabled. It is best, for starters, to replace all
  update() calls with set(). Recall, also, that direct assignment is now
  possible for adding new keywords to a header. So by and large the only
  reason to prefer using :meth:`Header.set` is its capability of inserting or
  moving a keyword to a specific location using the ``before`` or ``after``
  arguments.

- Slicing a Header with a slice index returns a new Header containing only
  those cards contained in the slice. As mentioned earlier, it used to be that
  slicing a Header returned a card list — something of a misfeature. In
  general, objects that support slicing ought to return an object of the same
  type when you slice them.

  Likewise, wildcard keywords used to return a CardList object — now they
  return a new Header similarly to a slice. For example::

      >>> header['NAXIS*']

  returns a new header containing only the NAXIS and NAXISn cards from the
  original header.

.. _duck typing: https://en.wikipedia.org/wiki/Duck_typing


Transition Tips
===============

The above may seem like a lot, but the majority of existing code using PyFITS
to manipulate headers should not need to be updated, at least not immediately.
The most common operations still work the same.

As mentioned above, it would be helpful to run your code with ``python -Wd`` to
enable deprecation warnings — that should be a good idea of where to look to
update your code.

If your code needs to be able to support older versions of PyFITS
simultaneously with PyFITS 3.1, things are slightly trickier, but not by
much — the deprecated interfaces will not be removed for several more versions
because of this.

- The first change worth making, which is supported by any PyFITS version in
  the last several years, is to remove any use of ``Header.has_key`` and
  replace it with ``keyword in header`` syntax. It is worth making this change
  for any dict as well, since `dict.has_key` is deprecated. Running the
  following regular expression over your code may help with most (but not all)
  cases::

      s/([^ ]+)\.has_key\(([^)]+)\)/\2 in \1/

- If possible, replace any calls to Header.update() with Header.set() (though
  do not bother with this if you need to support older PyFITS versions). Also,
  if you have any calls to Header.update() that can be replaced with simple
  subscript assignments (e.g., ``header['NAXIS'] = (2, 'Number of axes')``) do
  that too, if possible.

- Find any code that uses ``header.ascard`` or ``header.ascardlist()``. First
  ascertain whether that code really needs to work directly on Card objects.
  If that is definitely the case, go ahead and replace those with
  ``header.cards`` — that should work without too much fuss. If you do need to
  support older versions, you may keep using ``header.ascard`` for now.

- In the off chance that you have any code that slices a header, it is best to
  take the result of that and create a new Header object from it. For
  example::

      >>> new_header = fits.Header(old_header[2:])

  This avoids the problem that in PyFITS <= 3.0 slicing a Header returns a
  CardList by using the result to initialize a new Header object. This will
  work in both cases (in PyFITS 3.1, initializing a Header with an existing
  Header just copies it, à la `list`).

- As mentioned earlier, locate any code that deletes keywords with ``del`` and
  make sure they either look before they leap (``if keyword in header:``) or
  ask forgiveness (``try/except KeyError:``).

Other Gotchas
-------------

- As mentioned above, it is not necessary to enter ``print header`` to display
  a header in an interactive Python prompt. Entering ``>>> header``
  by itself is sufficient. Using ``print`` usually will *not* display the
  header readably, because it does not include line breaks between the header
  cards. The reason is that Python has two types of string representations.
  One is returned when a user calls ``str(header)``, which happens automatically
  when you ``print`` a variable. In the case of the Header class this actually
  returns the string value of the header as it is written literally in the
  FITS file, which includes no line breaks.

  The other type of string representation happens when one calls
  ``repr(header)``. The `repr` of an object is meant to be a useful
  string "representation" of the object; in this case the contents of the
  header but with line breaks between the cards and with the END card and
  trailing padding stripped off. This happens automatically when
  a user enters a variable at the Python prompt by itself without a ``print``
  call.

- The current version of the FITS Standard (3.0) states in section 4.2.1
  that trailing spaces in string values in headers are not significant and
  should be ignored. PyFITS < 3.1 *did* treat trailing spaces as significant.
  For example, if a header contained:

      KEYWORD1= 'Value    '

  then ``header['KEYWORD1']`` would return the string ``'Value    '`` exactly,
  with the trailing spaces intact. The new Header interface fixes this by
  automatically stripping trailing spaces, so that ``header['KEYWORD1']`` would
  return just ``'Value'``.

  There is, however, one convention used by the IRAF CCD mosaic task for
  representing its TNX World Coordinate System and ZPX World Coordinate System
  nonstandard WCS that uses a series of keywords in the form ``WATj_nnn``,
  which store a text description of coefficients for a nonlinear distortion
  projection. It uses its own microformat for listing the coefficients as a
  string, but the string is long, and thus broken up into several of these
  ``WATj_nnn`` keywords. Correct recombination of these keywords requires
  treating all whitespace literally. This convention either overlooked or
  predated the prescribed treatment of whitespace in the FITS standard.

  To get around this issue, a global variable ``fits.STRIP_HEADER_WHITESPACE``
  was introduced. Temporarily setting
  ``fits.STRIP_HEADER_WHITESPACE.set(False)`` before reading keywords affected
  by this issue will return their values with all trailing whitespace intact.

  A future version of PyFITS may be able to detect use of conventions like this
  contextually and behave according to the convention, but in most cases the
  default behavior of PyFITS is to behave according to the FITS Standard.
