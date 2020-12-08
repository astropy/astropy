.. include:: references.txt

.. _working_with_timezones:

Working with Timezones
**********************

It is sometimes needed to find the timezone of an observer from their
coordinates on Earth. Since timezones are defined legally for a site, they
often don't simply follow the longitude of that site, as such, in order
to figure out the timezone we need to consult a database.

For this purpose we are going to use the IANA time zone database, which is
part of Python's standard library with the `zoneinfo` module.

.. _prerequisites:

Prerequisites
=============
This functionality depends on the the `zoneinfo` module. This module is part
of Python >= 3.9 standard library. In case you are using a lower Python version,
you will need to install the `zoneinfo` module using the backports namespace
package.

Example
=======

..
  EXAMPLE START
  Determining the altitude/azimuth of a celestial object using the local observation time for a given geographical location

Start by importing the necessary modules:

.. code-block:: python

  >>> import datetime
  >>> from timezonefinder import TimezoneFinder
  >>> # in case we aren't on python 3.9 import zoneinfo from backports
  >>> try:
  >>>   from zoneinfo import ZoneInfo
  >>> except ImportError:
  >>>   from backports.zoneinfo import ZoneInf
  >>> import astropy.coordinates as coord
  >>> from astropy.time import Time
  >>> from astropy import units as 
  >>> import numpy as np


Define the target object:

.. code-block:: python

  >>> target = coord.SkyCoord.from_name('albireo')

Define the geographical location and get its zone info:
 
.. code-block:: python

  >>> loc = coord.EarthLocation.of_site('Apache Point Observatory')
  >>> tz_name = TimezoneFinder().timezone_at(lng=loc.lon.degree, lat=loc.lat.degree)
  >>> tz = ZoneInfo(tz_name)

Obtain the observation time from the local civil time and the geographical location:

.. code-block:: python

  >>> dt = datetime.datetime(2020, 5, 3, 18, 0, 0, tzinfo=tz)
  >>> obstime = Time(dt) + np.linspace(0, 12, 256) * u.hour

Extract observation data from the position of the object:

.. code-block:: python

  >>> altaz = coord.AltAz(location=loc, obstime=obstime)
  >>> airmass = target.transform_to(altaz).secz

..
  EXAMPLE END
