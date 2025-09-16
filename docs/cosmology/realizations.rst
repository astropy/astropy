.. _astropy-cosmology-realizations:

*******************************
Included Cosmology Realizations
*******************************

.. currentmodule:: astropy.cosmology.realizations

A number of preloaded cosmologies are available from analyses using the WMAP and
Planck satellite data. The full list of the predefined cosmologies is given by
``cosmology.realizations.available`` and summarized below:

===========  ============================== ====  ===== =======
Name         Source                         H0    Om    Flat
===========  ============================== ====  ===== =======
|WMAP1|      Spergel et al. 2003            72.0  0.257 Yes
|WMAP3|      Spergel et al. 2007            70.1  0.276 Yes
|WMAP5|      Komatsu et al. 2009            70.2  0.277 Yes
|WMAP7|      Komatsu et al. 2011            70.4  0.272 Yes
|WMAP9|      Hinshaw et al. 2013            69.3  0.287 Yes
|Planck13|   Planck Collab 2013, Paper XVI  67.8  0.307 Yes
|Planck15|   Planck Collab 2015, Paper XIII 67.7  0.307 Yes
|Planck18|   Planck Collab 2018, Paper VI   67.7  0.310 Yes
===========  ============================== ====  ===== =======

Currently, all are instances of |FlatLambdaCDM|. More details about exactly
where each set of parameters comes from are available in the docstring for each
object and described below::

  >>> from astropy.cosmology import WMAP7
  >>> print(WMAP7.__doc__)
  WMAP7 instance of FlatLambdaCDM cosmology
  (from Komatsu et al. 2011, ApJS, 192, 18, doi: 10.1088/0067-0049/192/2/18.
   Table 1 (WMAP + BAO + H0 ML).)


.. _astropy_cosmology_realizations_Planck18:

Planck 2018
===========

The parameters are from Planck Collaboration (2020) Table 2 (TT, TE, EE
+ lowE + lensing + BAO) [P18]_.

.. note::

  The Planck 2018 paper includes massive neutrinos in ``Om0`` but the Planck18
  object includes them in ``m_nu`` instead for consistency. Hence, the ``Om0``
  value in Planck18 differs slightly from the Planck 2018 paper but represents
  the same cosmological model.

Metadata
--------

==========  =====================================
Oc0         Omega cold dark matter at z=0
n           Density perturbation spectral index
sigma8      Density perturbation amplitude
tau         Ionisation optical depth
z_reion     Redshift of hydrogen reionisation
t0          Age of the universe in Gyr
reference   Reference for the parameters
==========  =====================================

References
----------
.. [P18] Planck Collaboration, et. al. (2020). Planck 2018 results. VI.
    Cosmological parameters. Astronomy \& Astrophysics, 641, A6.
    `<https://ui.adsabs.harvard.edu/abs/2020A%26A...641A...6P/abstract>`_


.. _astropy_cosmology_realizations_Planck15:

Planck 2015
===========

Parameters are from Planck Collaboration (2016) Paper XIII, Table 4 (TT, TE, EE + lowP + lensing + ext) [P15]_.

Metadata
--------

==========  =====================================
Oc0         Omega cold dark matter at z=0
n           Density perturbation spectral index
sigma8      Density perturbation amplitude
tau         Ionisation optical depth
z_reion     Redshift of hydrogen reionisation
t0          Age of the universe in Gyr
reference   Reference for the parameters
==========  =====================================

References
----------
.. [P15] Planck Collaboration, et. al. (2016). Planck 2015 results. XIII.
    Cosmological parameters. Astronomy \& Astrophysics, 594, A13.
    `<https://ui.adsabs.harvard.edu/abs/2016A%26A...594A..13P/abstract>`_


.. _astropy_cosmology_realizations_Planck13:

Planck 2013
===========

Parameters are from Planck Collaboration (2014) Paper XVI, Table 5 (Planck + WP
+ highL + BAO) [P13]_.

Metadata
--------

==========  =====================================
Oc0         Omega cold dark matter at z=0
n           Density perturbation spectral index
sigma8      Density perturbation amplitude
tau         Ionisation optical depth
z_reion     Redshift of hydrogen reionisation
t0          Age of the universe in Gyr
reference   Reference for the parameters
==========  =====================================

References
----------
.. [P13] Planck Collaboration, et. al. (2014). Planck 2013 results. XVI.
    Cosmological parameters. Astronomy \& Astrophysics, 571, A16.
    `<https://ui.adsabs.harvard.edu/abs/2014A%26A...571A..16P/abstract>`_


.. _astropy_cosmology_realizations_WMAP9:

WMAP 9 Year
===========

Parameters are from Hinshaw et al. (2013) Table 4 (WMAP9 + eCMB + BAO + H0) [WM9]_.

Metadata
--------

==========  =====================================
Oc0         Omega cold dark matter at z=0
n           Density perturbation spectral index
sigma8      Density perturbation amplitude
tau         Ionisation optical depth
z_reion     Redshift of hydrogen reionisation
t0          Age of the universe in Gyr
reference   Reference for the parameters
==========  =====================================

References
----------
.. [WM9] Hinshaw, G., Larson, D., Komatsu, E., Spergel, D., Bennett, C.,
    Dunkley, J., Nolta, M., Halpern, M., Hill, R., Odegard, N., Page, L., Smith,
    K., Weiland, J., Gold, B., Jarosik, N., Kogut, A., Limon, M., Meyer, S.,
    Tucker, G., Wollack, E., \& Wright, E. (2013). Nine-year Wilkinson Microwave
    Anisotropy Probe (WMAP) Observations: Cosmological Parameter Results. The
    Astrophysical Journal Supplement, 208(2), 19.
    `<https://ui.adsabs.harvard.edu/abs/2014A%26A...571A..16P/abstract>`_


.. _astropy_cosmology_realizations_WMAP7:

WMAP 7 Year
===========

Parameters are from Komatsu et al. (2011) Table 1 (WMAP + BAO + H0 ML) [WM7]_.

Metadata
--------

==========  =====================================
Oc0         Omega cold dark matter at z=0
n           Density perturbation spectral index
sigma8      Density perturbation amplitude
tau         Ionisation optical depth
z_reion     Redshift of hydrogen reionisation
t0          Age of the universe in Gyr
reference   Reference for the parameters
==========  =====================================

References
----------
.. [WM7] Komatsu, E., Smith, K., Dunkley, J., Bennett, C., Gold, B., Hinshaw,
    G., Jarosik, N., Larson, D., Nolta, M., Page, L., Spergel, D., Halpern, M.,
    Hill, R., Kogut, A., Limon, M., Meyer, S., Odegard, N., Tucker, G., Weiland,
    J., Wollack, E., \& Wright, E. (2011). Seven-year Wilkinson Microwave
    Anisotropy Probe (WMAP) Observations: Cosmological Interpretation. The
    Astrophysical Journal Supplement, 192(2), 18.
    `<https://ui.adsabs.harvard.edu/abs/2011ApJS..192...18K/abstract>`_


.. _astropy_cosmology_realizations_WMAP5:

WMAP 5 Year
===========

Parameters are from Komatsu et al. (2009) Table 1 (WMAP + BAO + SN ML) [WM5]_.

Metadata
--------

==========  =====================================
Oc0         Omega cold dark matter at z=0
n           Density perturbation spectral index
sigma8      Density perturbation amplitude
tau         Ionisation optical depth
z_reion     Redshift of hydrogen reionisation
t0          Age of the universe in Gyr
reference   Reference for the parameters
==========  =====================================

References
----------
.. [WM5] Komatsu, E., Dunkley, J., Nolta, M., Bennett, C., Gold, B., Hinshaw,
    G., Jarosik, N., Larson, D., Limon, M., Page, L., Spergel, D., Halpern, M.,
    Hill, R., Kogut, A., Meyer, S., Tucker, G., Weiland, J., Wollack, E., \&
    Wright, E. (2009). Five-Year Wilkinson Microwave Anisotropy Probe
    Observations: Cosmological Interpretation. The Astrophysical Journal
    Supplement, 180(2), 330-376.
    `<https://ui.adsabs.harvard.edu/abs/2009ApJS..180..330K/abstract>`_


.. _astropy_cosmology_realizations_WMAP3:

WMAP 3 Year
===========

Parameters are from Spergel et al. (2007) Table 6 (WMAP + SNGold) [WM3]_,
obtained from
https://lambda.gsfc.nasa.gov/product/map/dr2/params/lcdm_wmap_sngold.cfm.
``Tcmb0`` and ``Neff`` are the standard values as also used for |WMAP5|,
|WMAP7|, |WMAP9|.

.. note:: Pending WMAP team approval and subject to change.

Metadata
--------

==========  =====================================
Oc0         Omega cold dark matter at z=0
n           Density perturbation spectral index
sigma8      Density perturbation amplitude
tau         Ionisation optical depth
z_reion     Redshift of hydrogen reionisation
t0          Age of the universe in Gyr
reference   Reference for the parameters
==========  =====================================

References
----------
.. [WM3] Spergel, D., Bean, R., Doré, O., Nolta, M., Bennett, C., Dunkley, J.,
    Hinshaw, G., Jarosik, N., Komatsu, E., Page, L., Peiris, H., Verde, L.,
    Halpern, M., Hill, R., Kogut, A., Limon, M., Meyer, S., Odegard, N., Tucker,
    G., Weiland, J., Wollack, E., \& Wright, E. (2007). Three-Year Wilkinson
    Microwave Anisotropy Probe (WMAP) Observations: Implications for Cosmology.
    The Astrophysical Journal Supplement, 170(2), 377-408.
    `<https://ui.adsabs.harvard.edu/abs/2007ApJS..170..377S/abstract>`_


.. _astropy_cosmology_realizations_WMAP1:

WMAP 1 Year
===========

Parameters are from Spergel et al. (2003) Table 7 (WMAP + CBI + ACBAR + 2dFGRS
+ Lya) [WM1]_. ``Tcmb0`` and ``Neff`` are the standard values as also used for
|WMAP5|, |WMAP7|, |WMAP9|.

.. note:: Pending WMAP team approval and subject to change.

Metadata
--------

==========  =====================================
Oc0         Omega cold dark matter at z=0
n           Density perturbation spectral index
sigma8      Density perturbation amplitude
tau         Ionisation optical depth
z_reion     Redshift of hydrogen reionisation
t0          Age of the universe in Gyr
reference   Reference for the parameters
==========  =====================================

References
----------
.. [WM1] Spergel, D., Verde, L., Peiris, H., Komatsu, E., Nolta, M., Bennett,
    C., Halpern, M., Hinshaw, G., Jarosik, N., Kogut, A., Limon, M., Meyer, S.,
    Page, L., Tucker, G., Weiland, J., Wollack, E., \& Wright, E. (2003).
    First-Year Wilkinson Microwave Anisotropy Probe (WMAP) Observations:
    Determination of Cosmological Parameters. The Astrophysical Journal
    Supplement, 148(1), 175-194.
    `<https://ui.adsabs.harvard.edu/abs/2003ApJS..148..175S/abstract>`_


Reference/API
=============

.. py:currentmodule:: astropy.cosmology.realizations

.. automodapi:: astropy.cosmology.realizations
   :no-main-docstr:
   :include-all-objects:
   :noindex:
   :skip: default_cosmology
