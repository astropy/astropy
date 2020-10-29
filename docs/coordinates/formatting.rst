Formatting Coordinate Strings
*****************************

.. todo: @taldcroft should change this to start with a discussion of SkyCoord's capabilities

Getting a string representation of a coordinate is most powerfully
approached by treating the components (e.g., RA and Dec) separately.

Examples
--------

..
  EXAMPLE START
  Getting and Formatting String Representations of Coordinates

To get the string representation of a coordinate::

  >>> from astropy.coordinates import ICRS
  >>> from astropy import units as u
  >>> coo = ICRS(187.70592*u.degree, 12.39112*u.degree)
  >>> str(coo.ra) + ' ' + str(coo.dec)
  '187d42m21.312s 12d23m28.032s'

To get better control over the formatting, you can use the angles'
:meth:`~astropy.coordinates.Angle.to_string` method (see :doc:`angles` for
more). For example::

  >>> rahmsstr = coo.ra.to_string(u.hour)
  >>> str(rahmsstr)
  '12h30m49.4208s'
  >>> decdmsstr = coo.dec.to_string(u.degree, alwayssign=True)
  >>> str(decdmsstr)
  '+12d23m28.032s'
  >>> rahmsstr + ' ' + decdmsstr
  u'12h30m49.4208s +12d23m28.032s'

You can also use Python's `format` string method to create more complex
string expressions, such as IAU-style coordinates or even full sentences::

  >>> (f'SDSS J{coo.ra.to_string(unit=u.hourangle, sep="", precision=2, pad=True)}'
  ...  f'{coo.dec.to_string(sep="", precision=2, alwayssign=True, pad=True)}')
  'SDSS J123049.42+122328.03'
  >>> f'The galaxy M87, at an RA of {coo.ra.hour:.2f} hours and Dec of {coo.dec.deg:.1f} degrees, has an impressive jet.'
  'The galaxy M87, at an RA of 12.51 hours and Dec of 12.4 degrees, has an impressive jet.'

..
  EXAMPLE END
