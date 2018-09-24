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
 >>> from astropy.coordinates.tests.utils import randomly_sample_sphere
 >>> from astropy.time import Time
 >>> from astropy.units import u
 >>> import numpy as np

 >>> # 1000 random locations on the sky
 >>> ra, dec, _ = randomly_sample_sphere(1000)
 >>> coos = SkyCoord(ra, dec)

 >>> # 300 times over the space of 10 hours
 >>> times = Time.now() + np.linspace(-5, 5, 300)*u.hour

 >>> # note the use of broadcasting so that 300 times broadcast against 1000 positions
 >>> aa_frame = coord.AltAz(obstime=times[:, np.newaxis], location=EarthLocation.of_site('lapalma'))

 >>> # calculate alt-az of each object at each time.
 >>> aa_coos = coos.transform_to(aa_frame)