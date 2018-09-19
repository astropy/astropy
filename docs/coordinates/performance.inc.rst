.. note that if this is changed from the default approach of using an *include* 
   (in index.rst) to a separate performance page, the header needs to be changed
   from === to ***, the filename extension needs to be changed from .inc.rst to 
   .rst, and a link needs to be added in the subpackage toctree

.. _astropy-coordinates-performance:

Performance Tips
================

Here we provide some tips and tricks for how to optimize performance of code
using `astropy.coordinates`.

Use broadcasting to transform many SkyCoords into frames with vector properties
 >>> from astropy.coordinates import SkyCoord, EarthLocation
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> added import statement to correct coord error
 >>> from astropy import coordinates as coord
 >>> from astropy.coordinates.tests.utils import randomly_sample_sphere
 >>> from astropy.time import Time
<<<<<<< HEAD
 >>> from astropy import units as u
=======
 >>> from astropy.coordinates.tests.utils import randomly_sample_sphere
 >>> from astropy.time import Time
>>>>>>> Added broadcasting example to performance section of Coord docs
=======
 >>> from astropy.units import u
>>>>>>> Added units import line to fix pr errors
 >>> import numpy as np
 >>> # 1000 random locations on the sky
 >>> ra, dec, _ = randomly_sample_sphere(1000)
 >>> coos = SkyCoord(ra, dec)
 >>> # 300 times over the space of 10 hours
 >>> times = Time.now() + np.linspace(-5, 5, 300)*u.hour
 >>> # note the use of broadcasting so that 300 times broadcast against 1000 positions
<<<<<<< HEAD
<<<<<<< HEAD
 >>> lapalma = EarthLocation.from_geocentric(5327448.9957829, -1718665.73869569, 3051566.90295403, unit='m')
 >>> aa_frame = coord.AltAz(obstime=times[:, np.newaxis], location= lapalma)
=======
 >>> aa_frame = coord.AltAz(obstime=times[:, np.newaxis], location=EarthLocation.of_site('lapalma'))
>>>>>>> Added broadcasting example to performance section of Coord docs
=======
 >>> lapalma = EarthLocation.from_geocentric(5327448.9957829, -1718665.73869569, 3051566.90295403, unit='m')
 >>> aa_frame = coord.AltAz(obstime=times[:, np.newaxis], location= lapalma)
>>>>>>> added EarthLocation of site lapalma from geocentric
 >>> # calculate alt-az of each object at each time.
 >>> aa_coos = coos.transform_to(aa_frame)