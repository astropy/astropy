Formatting Coordinate Strings
-----------------------------

Getting a string representation of a coordinate is best approached by
treating the components (e.g., RA and Dec) separately.  For example::

  >>> from astropy.coordinates import ICRSCoordinates
  >>> c = ICRSCoordinates(187.70592, 12.39112, unit=(u.degree, u.degree))
  >>> str(c.ra) + ' ' + str(c.dec)
  '187d42m21.31200s 12d23m28.03200s'

To get better control over the formatting, you can use the angles' 
`~astropy.coordinates.angle.Angle.format` method (see :doc:`angles` for
more).  For example::

  >>> rahmsstr = c.ra.format(u.hour)
  >>>> rahmsstr
  '12h30m49.42080s'
  >>> decdmsstr = c.dec.format(u.degree,alwayssign=True)
  >>> decdmsstr
  '+12d23m28.03200s'
  >>> rahmsstr + ' ' + decdmsstr
  '12h30m49.42080s +12d23m28.03200s'
