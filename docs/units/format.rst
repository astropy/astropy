.. _astropy-units-format:

String Representations of Units and Quantities
**********************************************

Converting to Strings
=====================

You can control the way that |Quantity| and |Unit| objects are rendered as
strings using the `Python Format String Syntax
<https://docs.python.org/3/library/string.html#format-string-syntax>`_
(demonstrated below using `f-strings
<https://www.python.org/dev/peps/pep-0498/>`_).

For a |Quantity|, format specifiers that are names of `Built-In Formats`_ are
applied to the |Quantity| unit, and if possible also to the value. Format
specifiers for numerical values, like ``.3f``, will be applied to the
|Quantity| value, without affecting the unit. Finally, specifiers like
``^20s``, which would apply to a string, will be applied to the string
representation of the |Quantity| as a whole. Format specifiers that apply to
the unit part of a |Quantity| are also applicable to a |Unit| instance.

Examples
--------

.. EXAMPLE START: Converting Units to String Representations

To render |Quantity| or |Unit| objects as strings::

    >>> from astropy import units as u
    >>> q = 10.5 * u.km
    >>> q
    <Quantity  10.5 km>
    >>> f"{q}"
    '10.5 km'
    >>> f"{q:latex}"
    '$10.5 \\; \\mathrm{km}$'
    >>> f"{q:+.3f}"
    '+10.500 km'
    >>> f"{q:^20}"
    '        10.5         km'
    >>> f"{q:^20s}"
    '     10.5 km        '

To format both the value and the unit separately, you can access the |Quantity|
attributes within format strings::

    >>> q = 10.5 * u.km
    >>> q
    <Quantity  10.5 km>
    >>> f"{q.value:.3f} in {q.unit}"
    '10.500 in km'

This might not work well with LaTeX strings, in which case it would be better
to use the `Quantity.to_string() <astropy.units.Quantity.to_string()>`
method::

    >>> q = 1.2478e12 * u.pc/u.Myr
    >>> f"{q:latex}"  # Might not have the number of digits we would like
    '$1.2478 \\times 10^{12} \\; \\mathrm{\\frac{pc}{Myr}}$'
    >>> f"{q.value:.3e} {q.unit:latex}"  # The value is not in LaTeX
    '1.248e+12 $\\mathrm{\\frac{pc}{Myr}}$'
    >>> q.to_string(format="latex", precision=4)  # Right number of LaTeX digits
    '$1.248 \\times 10^{12} \\; \\mathrm{\\frac{pc}{Myr}}$'

Because |ndarray| does not accept most format specifiers, using specifiers like
``.3f`` will not work when applied to a |ndarray| or non-scalar |Quantity|. Use
:func:`numpy.array_str` instead. For instance::

    >>> import numpy as np
    >>> q = np.linspace(0,1,10) * u.m
    >>> f"{np.array_str(q.value, precision=1)} {q.unit}"  # doctest: +FLOAT_CMP
    '[0.  0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9 1. ] m'

Examine the NumPy documentation for more examples with :func:`numpy.array_str`.

.. EXAMPLE END

A |Unit|, or the unit part of a |Quantity|, can also be formatted in a number
of different styles. By default, the string format used is the "generic"
format, which is based on syntax of the `FITS standard
<https://fits.gsfc.nasa.gov/fits_standard.html>`_ format for representing
units, but supports all of the units defined within the :mod:`astropy.units`
framework, including user-defined units. The format specifier (and
`UnitBase.to_string() <astropy.units.core.UnitBase.to_string>`) functions also
take an optional parameter to select a different format::

    >>> q = 10 * u.km
    >>> f"{q:latex}"
    '$10 \\; \\mathrm{km}$'
    >>> fluxunit = u.erg / (u.cm ** 2 * u.s)
    >>> f"{fluxunit}"
    'erg / (s cm2)'
    >>> print(f"{fluxunit:unicode}")
    erg s⁻¹ cm⁻²
    >>> f"{fluxunit:latex}"
    '$\\mathrm{\\frac{erg}{s\\,cm^{2}}}$'
    >>> f"{fluxunit:>20s}"
    '       erg / (s cm2)'

The `UnitBase.to_string() <astropy.units.core.UnitBase.to_string>` method is an
alternative way to format units as strings, and is the underlying
implementation of the `format`-style usage::

    >>> fluxunit = u.erg / (u.cm ** 2 * u.s)
    >>> fluxunit.to_string('latex')
    '$\\mathrm{\\frac{erg}{s\\,cm^{2}}}$'

Converting from Strings
=======================

.. EXAMPLE START: Creating Units from Strings

Units can also be created from strings in a number of different
formats using the `~astropy.units.Unit` class::

  >>> u.Unit("m")
  Unit("m")
  >>> u.Unit("erg / (s cm2)")
  Unit("erg / (s cm2)")
  >>> u.Unit("erg.s-1.cm-2", format="cds")
  Unit("erg / (s cm2)")

It is also possible to create a scalar |Quantity| from a string::

    >>> u.Quantity("3m/s")
    <Quantity 3. m / s>

.. note::

   Converting from strings requires the use of a specialized parser for the
   unit language, which results in a performance penalty. It is much faster to
   use |Unit| objects directly (e.g., ``unit = u.degree / u.minute``) instead
   of via string parsing (``unit = u.Unit('deg/min')``). This parser is very
   useful, however, if your unit definitions are coming from a file format such
   as FITS or VOTable.

.. EXAMPLE END

Built-In Formats
================

`astropy.units` includes support for parsing and writing the following
formats:

  - ``"fits"``: This is the format defined in the Units section of the
    `FITS Standard <https://fits.gsfc.nasa.gov/fits_standard.html>`__.
    Unlike the "generic" string format, this will only accept or
    generate units defined in the FITS standard.

  - ``"vounit"``: The `Units in the VO 1.0
    <http://www.ivoa.net/documents/VOUnits/>`__ standard for
    representing units in the VO. Again, based on the FITS syntax,
    but the collection of supported units is different.

  - ``"cds"``: `Standards for astronomical catalogues from Centre de
    Données astronomiques de Strasbourg
    <https://vizier.unistra.fr/vizier/doc/catstd-3.2.htx>`_: This is the
    standard used by `Vizier tables <https://vizier.unistra.fr/>`__,
    as well as what is used by VOTable versions 1.3 and earlier.

  - ``"ogip"``: A standard for storing units as recommended by the
    `Office of Guest Investigator Programs (OGIP)
    <https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`_.

`astropy.units` is also able to write, but not read, units in the
following formats:

  - ``"latex"``: Writes units out using LaTeX math syntax using the
    `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_
    recommendations for unit presentation. This format is
    automatically used when printing a unit in the |IPython| notebook::

        >>> f"{fluxunit:latex}"
        '$\\mathrm{\\frac{erg}{s\\,cm^{2}}}$'

    which renders as

    .. math::

       \mathrm{\frac{erg}{s\,cm^{2}}}

  - ``"latex_inline"``: Writes units out using LaTeX math syntax using the
    `IAU Style Manual
    <https://www.iau.org/static/publications/stylemanual1989.pdf>`_
    recommendations for unit presentation, using negative powers instead of
    fractions, as required by some journals (e.g., `Apj and AJ
    <https://journals.aas.org/manuscript-preparation/>`_).
    Best suited for unit representation inline with text::

        >>> fluxunit.to_string('latex_inline')
        '$\\mathrm{erg\\,s^{-1}\\,cm^{-2}}$'

    which renders as

    .. math::

       \mathrm{erg\,s^{-1}\,cm^{-2}}

  - ``"console"``: Writes a representation of the unit useful for
    display in a text console::

      >>> print(fluxunit.to_string('console'))
       erg s^-1 cm^-2

    It is also possible to use a fraction, either on a single line,

      >>> print(fluxunit.to_string('console', fraction='inline'))
      erg / (s cm^2)

    or using a multiline representation:

      >>> print(fluxunit.to_string('console', fraction='multiline'))
       erg
      ------
      s cm^2

  - ``"unicode"``: Same as ``"console"``, except uses Unicode
    characters::

      >>> print(u.Ry.decompose().to_string('unicode'))  # doctest: +FLOAT_CMP
      2.1798724×10⁻¹⁸ m² kg s⁻²
      >>> print(u.Ry.decompose().to_string('unicode', fraction=True))  # doctest: +FLOAT_CMP
      2.1798724×10⁻¹⁸ m² kg / s²
      >>> print(u.Ry.decompose().to_string('unicode', fraction='multiline'))  # doctest: +FLOAT_CMP
                      m² kg
      2.1798724×10⁻¹⁸ ─────
                       s²

.. _astropy-units-format-unrecognized:

Dealing with Unrecognized Units
===============================

Since many files found in the wild have unit strings that do not
correspond to any given standard, `astropy.units` also has a
consistent way to store and pass around unit strings that did not
parse.  In addition, it provides tools for transforming non-standard,
legacy or misspelt unit strings into their standardized form,
preventing the further propagation of these unit strings.

By default, passing an unrecognized unit string raises an exception::

  >>> # The FITS standard uses 'angstrom', not 'Angstroem'
  >>> u.Unit("Angstroem", format="fits")
  Traceback (most recent call last):
    ...
  ValueError: 'Angstroem' did not parse as fits unit: At col 0, Unit
  'Angstroem' not supported by the FITS standard. Did you mean Angstrom
  or angstrom? If this is meant to be a custom unit, define it with
  'u.def_unit'. To have it recognized inside a file reader or other
  code, enable it with 'u.add_enabled_units'. For details, see
  https://docs.astropy.org/en/latest/units/combining_and_defining.html

However, the `~astropy.units.Unit` constructor has the keyword
argument ``parse_strict`` that can take one of three values to control
this behavior:

  - ``'raise'``: (default) raise a :class:`ValueError`.

  - ``'warn'``: emit a :class:`~astropy.units.UnitsWarning`, and return an
    `~astropy.units.UnrecognizedUnit` instance.

  - ``'silent'``: return an `~astropy.units.UnrecognizedUnit`
    instance.

By either adding additional unit aliases for the misspelt units with
:func:`~astropy.units.set_enabled_aliases` (e.g., 'Angstroms' for 'Angstrom';
as demonstrated below), or defining new units via
:func:`~astropy.units.def_unit` and :func:`~astropy.units.add_enabled_units`,
we can use ``parse_strict='raise'`` to rapidly find issues with the units used,
while also being able to read in older datasets where the unit usage may have
been less standard.


Examples
--------

.. EXAMPLE START: Define Aliases for Units

To set unit aliases, pass :func:`~astropy.units.set_enabled_aliases` a
:class:`dict` mapping the misspelt string to an astropy unit. The following
code snippet shows how to set up Angstroem -> Angstrom::

    >>> u.set_enabled_aliases({"Angstroem": u.Angstrom})
    <astropy.units.core._UnitContext object at 0x...>
    >>> u.Unit("Angstroem")
    Unit("Angstrom")
    >>> u.Unit("Angstroem") == u.Angstrom
    True

You can also set multiple aliases up at once or add to existing ones::

    >>> u.set_enabled_aliases({"Angstroem": u.Angstrom, "Angstroms": u.Angstrom})
    <astropy.units.core._UnitContext object at 0x...>
    >>> u.add_enabled_aliases({"angstroem": u.Angstrom})
    <astropy.units.core._UnitContext object at 0x...>
    >>> u.Unit("Angstroem") == u.Unit("Angstroms") == u.Unit("angstroem") == u.Angstrom
    True

The aliases can be reset by passing an empty dictionary::

    >>> u.set_enabled_aliases({})
    <astropy.units.core._UnitContext object at 0x...>

You can use both :func:`~astropy.units.set_enabled_aliases` and
:func:`~astropy.units.add_enabled_aliases` as a `context manager
<https://docs.python.org/3/reference/datamodel.html#context-managers>`_,
limiting where a particular alias is used::

    >>> with u.add_enabled_aliases({"Angstroem": u.Angstrom}):
    ...     print(u.Unit("Angstroem") == u.Angstrom)
    True
    >>> u.Unit("Angstroem") == u.Angstrom
    Traceback (most recent call last):
      ...
    ValueError: 'Angstroem' did not parse as unit: At col 0, Angstroem is not
    a valid unit. Did you mean Angstrom, angstrom, mAngstrom or mangstrom? If
    this is meant to be a custom unit, define it with 'u.def_unit'. To have it
    recognized inside a file reader or other code, enable it with
    'u.add_enabled_units'. For details, see
    https://docs.astropy.org/en/latest/units/combining_and_defining.html

.. EXAMPLE END

.. EXAMPLE START: Using `~astropy.units.UnrecognizedUnit`

To pass an unrecognized unit string::

   >>> x = u.Unit("Angstroem", format="fits", parse_strict="warn")  # doctest: +SHOW_WARNINGS
   UnitsWarning: 'Angstroem' did not parse as fits unit: At col 0, Unit
   'Angstroem' not supported by the FITS standard. Did you mean Angstrom or
   angstrom? If this is meant to be a custom unit, define it with 'u.def_unit'.
   To have it recognized inside a file reader or other code, enable it with
   'u.add_enabled_units'. For details, see
   https://docs.astropy.org/en/latest/units/combining_and_defining.html

This `~astropy.units.UnrecognizedUnit` object remembers the
original string it was created with, so it can be written back out,
but any meaningful operations on it, such as converting to another
unit or composing with other units, will fail.

   >>> x.to_string()
   'Angstroem'
   >>> x.to(u.km)
   Traceback (most recent call last):
     ...
   ValueError: The unit 'Angstroem' is unrecognized.  It can not be
   converted to other units.
   >>> x / u.m
   Traceback (most recent call last):
     ...
   ValueError: The unit 'Angstroem' is unrecognized, so all arithmetic
   operations with it are invalid.

.. EXAMPLE END
