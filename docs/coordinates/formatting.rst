Formatting Coordinate Strings
-----------------------------

.. todo: @taldcroft should change this to start with a discussion of SkyCoord's capabilities

Getting a string representation of a coordinate is most powerfuly
approached by treating the components (e.g., RA and Dec) separately.
For example::

  >>> from astropy.coordinates import ICRS
  >>> from astropy import units as u
  >>> c = ICRS(187.70592*u.degree, 12.39112*u.degree)
  >>> str(c.ra) + ' ' + str(c.dec)
  '187d42m21.312s 12d23m28.032s'

To get better control over the formatting, you can use the angles'
:meth:`~astropy.coordinates.Angle.to_string` method (see :doc:`angles` for
more).  For example::

  >>> rahmsstr = c.ra.to_string(u.hour)
  >>> str(rahmsstr)
  '12h30m49.4208s'
  >>> decdmsstr = c.dec.to_string(u.degree, alwayssign=True)
  >>> str(decdmsstr)
  '+12d23m28.032s'
  >>> rahmsstr + ' ' + decdmsstr
  u'12h30m49.4208s +12d23m28.032s'

You can also use python's `format` string method to create more complex
string expressions, such as IAU-style coordinates or even full sentences::

  >>> 'SDSS J{0}{1}'.format(c.ra.to_string(sep='', precision=2, pad=True), c.dec.to_string(sep='', precision=2, alwayssign=True, pad=True))
  'SDSS J1874221.31+122328.03'
  >>> 'The galaxy M87, at an RA of {0.ra.deg:.1f} and Dec of {0.dec.deg:.1f} degrees, has an impressive jet.'.format(c)
  'The galaxy M87, at an RA of 187.7 and Dec of 12.4 degrees, has an impressive jet.'
