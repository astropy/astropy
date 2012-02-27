# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .. import cosmology

def test_cosmology():
    cosmo = cosmology.Cosmology(H0=70, Om=0.27 ,Ol=0.73)
    z = 1
    # Test values were taken from the following web cosmology
    # calculators on 27th Feb 2012:

    # Wright: http://www.astro.ucla.edu/~wright/CosmoCalc.html
    # Kempner: http://www.kempner.net/cosmic.php
    # iCosmos: http://www.icosmos.co.uk/index.html

    # The order of values below is Wright, Kempner, iCosmos'
    assert np.allclose(cosmo.dc(z), [3364.5, 3364.8, 3364.7988], rtol=1e-4)
    assert np.allclose(cosmo.da(z), [1682.3, 1682.4, 1682.3994], rtol=1e-4)
    assert np.allclose(cosmo.dl(z), [6729.2, 6729.6, 6729.5976], rtol=1e-4)
    assert np.allclose(cosmo.tl(z), [7.841, 7.84178, 7.843],  rtol=1e-3)

