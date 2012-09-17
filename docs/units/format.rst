Unit formats
============

Units can be created from strings using the `Unit` class::

  >>> u.Unit("m")
  Unit("m")
  >>> u.Unit("erg / (s cm2)")
  Unit("erg / (s cm2)")

Units can be converted to strings using the `to_string` method::

  >>> fluxunit = u.erg / (u.cm ** 2 * u.s)
  >>> fluxunit.to_string()
  u'erg / (s cm2)'

By default, the string format used is referred to as the "generic"
format, which is based on syntax of the FITS standard's format for
representing units, but supports all of the units defined within the
`astropy.units` framework, including user-defined units.  The `Unit`
and `to_string` functions also take an optional `format` parameter to
select a different format.  This parameter may be either a string or a
`astropy.units.format.Base` instance.

Built-in formats
----------------

`astropy.units` includes support for parsing and writing the following
formats:

  - ``"fits"``: This is the format defined in the Units section of the
    `FITS Standard <http://fits.gsfc.nasa.gov/fits_standard.html>`_.
    Unlike the "generic" string format, this will only accept or
    generate units defined in the FITS standard.

  - ``"vounit"``: The `proposed IVOA standard
    <http://www.ivoa.net/Documents/VOUnits/>`_ for representing units
    in the VO.  Again, based on the FITS syntax, but the collection of
    supported units is different.

.. These are to-be-implemented

  - OGIP Units: A standard for storing units in `OGIP FITS files
    <http://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/general/ogip_93_001/>`_.

  - `Standards for astronomical catalogues
    <http://cds.u-strasbg.fr/doc/catstd-3.2.htx>`_: This is the
    standard used, for example, by VOTable versions 1.2 and earlier.

`astropy.units` also is able to write units in the following formats:

  - ``"latex"``: Writes units out using LaTeX math syntax using the
    `IAU Style Manual
    <http://www.iau.org/static/publications/stylemanual1989.pdf>`_
    recommendations for unit presentation.  This format is
    automatically used when printing a unit in the IPython notebook::

      >>> fluxunit

    .. math::

       \mathrm{\frac{erg}{s\ cm^{2}}}

  - ``"console"``: Writes a multi-line representation of the unit
    useful for display in a text console::

      >>> print fluxunit.to_string('console')
       erg
      ------
      s cm^2

  - ``"unicode"``: Same as ``"console"``, except uses Unicode
    characters::

      >>> print u.Ry.decompose().to_string('unicode')
                 m² kg
      2.18×10-¹⁸ ─────
                  s²
